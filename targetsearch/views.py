from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse, Http404
from django.urls import reverse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from tools.models import *
from django.contrib import messages
from collections import OrderedDict
from django.template import Context, Engine
from guest.decorators import guest_allowed
import sys
import traceback

from .helpers import *
from . import ts_graph

from django.conf import settings
from django.views.decorators.cache import cache_page

import os
import csv
import json

def readSources(type):
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,type+"_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[key.strip()] = val.strip()
    return sources

def detailPage(request,id):
    targetInfo=None
    compoundInfo=None
    id_type = None
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    activity_info= None
    activity_matches= None


    #print("rendering detail page for "+str(id))

    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    else:
        id_type = getEntityType(id)

    if id_type == None or id_type == "unknown":
        return render(request,'detail-select-database.html',{
            'given_id': id,
            'compoundDbs' : compoundDbs,
            'proteinDbs':proteinDbs,
        })


    try:
        if id_type == "compound" and 'compound_source_id' in request.GET:
            source_id = request.GET['compound_source_id']
            #get mapping dict, fetch keys, then get either first element, or None
            id=next(iter(mapToChembl([id], source_id).keys()), None)
            #print("new compound id: "+str(id)+", source_id: "+str(source_id))
        elif id_type == "target" and 'target_source_id' in request.GET:
            source_id = request.GET['target_source_id']
            #get mapping dict, fetch keys, then get either first element, or None
            id=next(iter(mapToUniprot([id], source_id).keys()), None)
            #print("new target id: "+str(id)+", source_id: "+str(source_id))

        if id == None:
            raise Exception("ID not found in selected database")
    except Exception as e:
        message = "Sorry, I could not convert that ID to a ChEMBL compound or target. Error was: "+str(e)
        return render(request,'detail-select-database.html',{
            'given_id': id,
            'compoundDbs' : compoundDbs,
            'proteinDbs':proteinDbs,
            'message': message,
        })


    myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, [id])
    annotation_info = myAnnotationSearch.table_info
    annotation_matches = myAnnotationSearch.get_grouped_results()


    drugind_json = drugIndicationData(myAnnotationSearch.drugind_objs);

    if id_type == "target":
        targetInfo = targetSummaryInfo(id)

    if id_type == "compound":
        compoundInfo = compoundSummaryInfo(id);
        # activity search is often very slow, so disable for now
        #myActivitySearch = ActivitySearch(id_type, [id])
        #activity_info = myActivitySearch.table_info
        #activity_matches = myActivitySearch.get_grouped_results()


    return render(request,'detail.html', {
        'given_id':id,
        'type': id_type,
        'annotation_info' : annotation_info,
        'annotation_matches' : annotation_matches,
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'drugind_json' : json.dumps(drugind_json),
        'target_info': targetInfo,
        'compound_info': compoundInfo,
        })

@guest_allowed
def newTS(request,initial_ids=""):
    # Default local variables
    query_submit = False
    message = None
    annotation_info = None
    annotation_matches = None
    drugind_json = None
    activity_info = None
    activity_matches = None
    homolog_type_value= None
    allTags = Tag.allUserTagNames(request.user)
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    defaultCompoundDb = "1"
    defaultProteinDb = "ACC+ID"
    groupingCol = 0
    similarityJobs = [(job.id,str(job)) for job in Job.objects.filter(user=request.user,application__category_id=5)]

    # Default GET request variables
    id_type = 'compound'
    ids = list()
    include_activity = False
    source_id = 1
    similarity_job_id=-1

    # Retrieve GET request variables
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    #see if we can identify the entity type of any initial_ids
    elif initial_ids != "": 
        ii_array = initial_ids.split(",") 
        if len(ii_array) > 0:
            type_guess = getEntityType(ii_array[0])
            if type_guess != "unknown":
                id_type = type_guess


    if 'ids' in request.GET:
        if id_type == 'homolog-target':
            ids = request.GET.getlist('ids')
        else:
            ids = list(request.GET['ids'].split())
        initial_ids = " ".join(ids)
    if 'include_activity' in request.GET:
        include_activity = True
    if 'tags' in request.GET:
        for c in [compound.cid for compound in Compound.byTagNames(request.GET.getlist("tags"),request.user)]:
            ids.append(c)
    if 'source_id' in request.GET:
        source_id = request.GET['source_id']
    if 'similarity_job_id' in request.GET:
        similarity_job_id = int(request.GET['similarity_job_id'])

    # Generate content
    try:
        idMapping={}
        if id_type == 'compound' and source_id != '1':
            idMapping= mapToChembl(ids, source_id)
            ids = list(idMapping.keys())
        elif id_type == 'target' and source_id != 'ACC':
            idMapping= mapToUniprot(ids, source_id)
            ids = list(idMapping.keys())
        elif id_type == 'homolog-target':
            #full_ids = ids # Full context for homolog handling
            #ids = [ i.split(',')[2] for i in full_ids ] # Bare Accession IDs for old target-search

            # Context dictionary for homolog handling
            # ex:
            # homolog_context = {
            #     'P29274': {'paralog': 'P30542'},
            #     'P29275': {'paralog': 'P30542'},
            #     'P0DMS8': {'paralog': 'P30542'},
            #     }
            homolog_context = dict()
            for i in ids:
                [relation, src_id, homolog_id] = i.split(',')
                if homolog_id not in homolog_context.keys():
                    homolog_context[homolog_id] = dict()
                if src_id == homolog_id:
                    continue # Don't bother with "X is a homolog of X"
                if relation not in homolog_context[homolog_id].keys():
                    homolog_context[homolog_id][relation] = set()
                homolog_context[homolog_id][relation].add(src_id)

            ids = list(homolog_context.keys()) # Bare Accession IDs for old target-search

            # Prepare homolog relation descriptions
            homolog_desc = dict()
            for homolog_id, relations in homolog_context.items():
                desc_parts = list()
                for relation, src_ids in sorted(relations.items()):
                    desc_parts.append("{} of {}".format(relation, ', '.join(sorted(list(src_ids)))))
                homolog_desc[homolog_id] = '; '.join(desc_parts)

        if len(ids) != 0:
            query_submit = True
            queryIdCol = {
                    "id": "query_id",
                    "sql": None,
                    "table": "Query ID",
                    "name": "Query ID",
                    "desc": "Original compound ID prior to ChEMBL conversion",
                    "visible": True,
                    }
            headerTextCol = {
                "id": "header_text",
                "sql": None,
                "table": "Header Text",
                "name": "Header Text",
                "desc": "Description text to show in row-group headers (i.e. source query IDs, translations, etc.)",
                "visible": False,
                }
            originalQueryCol = {
                    "id": "original_query_id",
                    "sql": None,
                    "table": "Original Query ID",
                    "name": "Original Query ID",
                    "desc": "The compound that the current query compound was originally derived from, based on similarity",
                    "visible": True,
                    }

            myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, ids)
            annotation_info = myAnnotationSearch.table_info
            annotation_matches = myAnnotationSearch.get_grouped_results()

            drugind_json = drugIndicationData(myAnnotationSearch.drugind_objs);

            # Exclude ActivitySearch from search-by-target by default
            if id_type in ['target', 'homolog-target'] and not include_activity:
                activity_info = None
                activity_matches = None
            else:
                myActivitySearch = ActivitySearch(id_type, ids)
                activity_info = myActivitySearch.table_info
                activity_matches = myActivitySearch.get_grouped_results()

            if len(idMapping) != 0:
                groupingCol += 1
                addMappedQueryColumn(idMapping,queryIdCol,
                        annotation_info,annotation_matches,
                        activity_info,activity_matches)

            if similarity_job_id != -1:
                similarityMapping = readSimilarityMappingData(request.user,similarity_job_id)
                if len(idMapping) != 0: # need to compose our mapping with previous mapping
                    similarityMapping = composeMaps(idMapping,similarityMapping)
                #print("similarity mapping: \n"+str(similarityMapping))
                if len(similarityMapping) != 0:
                    groupingCol += 1
                    addMappedQueryColumn(similarityMapping, originalQueryCol,
                            annotation_info,annotation_matches,
                            activity_info,activity_matches)

            #if similarity_job_id != -1, then read job_<similarity_job_id>
            # map chembl id (2nd column) back to the original compound id (first column)
            # insert new column to show original compound id.

        if ts_paralog_cache():
            homolog_type_value = 'paralog-cache'
        else:
            homolog_type_value = 'paralog'

        initial_ids = initial_ids.replace(","," ")

    except Exception as e:
        print("exception in newTS:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        message = str(e)

    ids_resubmit = '+'.join(ids)

    context = {
        'query_submit' : query_submit,
        'message' : message,
        'id_type' : id_type,
        'annotation_info' : annotation_info,
        'annotation_matches' : annotation_matches,
        'drugind_json' : json.dumps(drugind_json),
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'tags' : allTags,
        'compoundDbs' : compoundDbs,
        'defaultCompoundDb': defaultCompoundDb,
        'proteinDbs': proteinDbs,
        'defaultProteinDb': defaultProteinDb,
        'groupingCol' : groupingCol,
        'similarityJobs': similarityJobs,
        'homolog_type_value': homolog_type_value,
        'initial_ids': initial_ids,
        'ids_resubmit' : ids_resubmit,
        }

    return render(request, 'targetsearch/new_ts.html', context)

def addMappedQueryColumn(mapping,columnDefinition,annotation_info, annotation_matches,activity_info,activity_matches):

    def addQueryCol(matches):
        for key, rowGroup in matches.items():
            if key in mapping:
                queryKey = mapping[key]
                for i in range(len(rowGroup)):
                    matches[key][i] = OrderedDict([(columnDefinition["id"],queryKey)] + [ item for item in rowGroup[i].items()] )


    annotation_info.insert(0,columnDefinition)

    addQueryCol(annotation_matches)
    if activity_matches != None:
        activity_info.insert(0,columnDefinition)
        addQueryCol(activity_matches)


def readSimilarityMappingData(user,job_id):
    job = Job.objects.get(id=job_id, user=user)
    mapping = {}

    try:
        f = open(job.output, 'r')
        csvinput = csv.reader(f, delimiter=' ')
        for line in csvinput:
            #map from target back to query
            # multiple queries may map to the same target
            queries = mapping.setdefault(line[1],[]);
            queries.append(line[0]);
    finally:
        f.close()

    for key, value in mapping.items():
        mapping[key] = ",".join(value)

    return mapping


## map1: x->y and map2: y -> z, return x->z
def composeMaps(map1,map2):
    composedMap = {}
    for key1,value1 in map1.items():
        if value1 in map2:
            composedMap[key1] = map2[value1]
    return composedMap

def tableHtml(search_obj, header=True, footer=False, table_class=""):
    """Takes a SearchBase object and returns an HTML table"""

    error = None
    table_info = None
    table_data = None

    engine = Engine.get_default()
    template = engine.get_template('targetsearch/table.html')

    try:
        table_info = search_obj.table_info
        table_data = search_obj.get_results()
    except Exception as e:
        error = str(e)

    context = {
        'error' : error,
        'header' : header,
        'footer' : footer,
        'table_info' : table_info,
        'table_data' : table_data,
        'table_class' : table_class,
    }

    return template.render(Context(context))

def drugIndTable(request):
    if 'molregno' in request.GET:
        molregno = request.GET['molregno']
    else:
        error = "drugIndTable: Missing molregno in GET variables"
        return render(request, 'targetsearch/table.html', {'error':error})

    myDrugIndSearch = DrugIndicationSearch(molregno)
    html = tableHtml(myDrugIndSearch)

    return HttpResponse(html)

def compoundNames(request,query):
    names = compoundNameAutocomplete(query)
    data = [ {'chembl_id': n.chembl_id, 'name': n.synonyms} for n in names]
    return JsonResponse(data,safe=False)

def targetNames(request, query):
    names = targetNameAutocomplete(query)
    data = [ {'accession_id': n.accession,
              'name': (n.description+' ('+n.organism+')')} for n in names ]
    return JsonResponse(data, safe=False)

def drugIndicationData(drugind_objs):
    # Generate Drug Indication JSON data
    drugind_json = dict()
    for chembl_id, drugind_obj in drugind_objs.items():
        colnames = [ {'title': c['name']} for c in drugind_obj.table_info ]

        data = list()
        for row in drugind_obj.get_results():
            data.append(list(row.values()))

        drugind_json[chembl_id] = dict()
        drugind_json[chembl_id]['colnames'] = colnames
        drugind_json[chembl_id]['data'] = data
    return drugind_json

@cache_page(60 * 120)
def chemblPNG(request, chembl_id):
    try:
        img = getChemblPNG(chembl_id, mwt_limit=2000, shrink=False)
        return HttpResponse(img, content_type='image/png')
    except Exception as e:
        raise Http404(str(e))

@cache_page(60 * 120)
def chemblSVG(request, chembl_id):
    try:
        img = getChemblSVG(chembl_id)
        return HttpResponse(img, content_type='image/svg+xml')
    except Exception as e:
        raise Http404(str(e))

def ajax(request, action):
    def die(msg):
        ajaxResponse = { 'success' : False, 'message' : msg }
        return JsonResponse(ajaxResponse)

    if action == 'homologs':
        try:
            ids = request.POST['ids']
            homolog_type = request.POST['homolog_type']
            result = homologsAjax(ids.split(), homolog_type)
            return JsonResponse(result)
        except Exception as e:
            return die(str(e))
    else:
        return die('Unknown action: {}'.format(action))

def extAnnoByChembl(request, chembl_id, db=None):
    check_avail = checkExtAnno(chembl_id)
    if db is None:
        data = None
    else:
        data = getExtAnno(chembl_id, db)

    context = {
        'check_avail': check_avail,
        'data': data,
        'chembl_id': chembl_id,
        'db': db,
    }

    return render(request, 'targetsearch/extannobychembl.html', context)

def tsGraph(request, id_type, table_name, ids):
    id_list = ids.split('+')

    if table_name == 'annotation':
        anno_data = ts_graph.get_anno_data(id_type, id_list)
        (nodes_data, edges_data) = ts_graph.make_anno_graph(anno_data)
    elif table_name == 'activity':
        act_data = ts_graph.get_act_data(id_type, id_list)
        (nodes_data, edges_data) = ts_graph.make_act_graph(act_data)

    context = {
        'id_type' : id_type,
        'table_name' : table_name,
        'ids' : ids,
        'ids_type' : str(type(ids)),
        'id_list' : id_list,
        'nodes_data' : nodes_data,
        'edges_data' : edges_data,
    }

    return render(request, 'targetsearch/ts_graph.html', context)

def tsFilter_debug(ids):
    myActivitySearch = ActivitySearch('compound', ids)
    activity_info = myActivitySearch.table_info
    activity_matches = myActivitySearch.get_results()

    #for m in activity_matches:
    #    m.pop('activity__component_sequences__sequence', None)

    cols = list()
    for k in activity_matches[0].keys():
        cols.append(k)

    #stdtypes_set = set()
    #for m in activity_matches:
    #    stdtypes_set.add(m['activity__activities__standard_type'])
    #stdtypes = list(stdtypes_set)
    #stdtypes.sort()

    # Standard Types can have different units, so we must use a (type, unit)
    # tuple to distinguish. Ugh.
    #stdtype_tuples_set = set()
    #for m in activity_matches:
    #    t = m['activity__activities__standard_type']
    #    u = m['activity__activities__standard_units']
    #    stdtype_tuples_set.add( (t,u) )
    #stdtype_tuples = list(stdtype_tuples_set)
    #stdtype_tuples.sort()

    # Standard Types can have different units, so we must use a (type, unit)
    # tuple to distinguish. Ugh.
    stdtype_tuples_set = set()
    stdval_low = dict()     # Lowest std value for a given (type, unit) tuple
    stdval_high = dict()    # Highest std value for a given (type, unit) tuple
    stdval_null = dict()    # Presence of NULL for a given (type, unit) tuple

    # Go thru every match to find all (t,u) tuples and lowest/highest values
    for m in activity_matches:
        t = m['activity__activities__standard_type']
        u = m['activity__activities__standard_units'] or 'NULL'
        v = m['activity__activities__standard_value']
        tu = (t,u)

        stdtype_tuples_set.add(tu)

        if v is None:
            stdval_null[tu] = True
        else:
            if tu not in stdval_low:
                stdval_low[tu] = v
            elif v < stdval_low[tu]:
                stdval_low[tu] = v
            if tu not in stdval_high:
                stdval_high[tu] = v
            elif v > stdval_high[tu]:
                stdval_high[tu] = v

    #stdtype_tuples = list()
    #for tu in stdtype_tuples_set:
    #    if tu in stdval_low and tu in stdval_high:
    #        stdtype_tuples.append(tu)
    #    else:
    #        t, u = tu
    #        stdtype_tuples.append( (t, 'YEET') )

    stdtype_tuples = list(stdtype_tuples_set)
    stdtype_tuples.sort()

    cols_min = [
        'activity__chembl_id_lookup__chembl_id',
        'activity__chembl_id_lookup__entity_type',
        'activity__activities__standard_type',
        'activity__activities__standard_relation',
        'activity__activities__standard_value',
        'activity__activities__standard_units',
        'activity__assays__assay_id',
        'activity__assays__doc_id',
        'activity__assays__description',
        'activity__assays__chembl_id',
        'activity__component_sequences__accession'
        ]

    context = {
        'ids' : ids,
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'cols' : cols,
        'cols_min' : cols_min,
        'stdtype_tuples' : stdtype_tuples,
        'stdval_low' : stdval_low,
        'stdval_high' : stdval_high,
        'stdval_null' : stdval_null,
    }

    return context

def tsFilter(request):
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    #else:
    #    ids = None

    context = tsFilter_debug(ids)

    return render(request, 'targetsearch/filter.html', context)

def tsFilter2_debug(ids):
    # Get annotation data, just as before
    myAnnoSearch = AnnotationWithDrugIndSearch('compound', ids)
    anno_info = myAnnoSearch.table_info
    anno_matches = myAnnoSearch.get_results()

    # Remove show_extanno. We don't need it here, plus it can't be JSONed
    for m in anno_matches:
        m.pop('annotation__extanno__show_extanno')

    cols = list()
    for k in anno_matches[0].keys():
        cols.append(k)

    # Minimal set of columns for a compact table
    cols_basic = [
        'annotation__chembl_id_lookup__chembl_id',
        'annotation__molecule_dictionary__pref_name',
        'annotation__drug_mechanism__mechanism_of_action',
        'annotation__component_sequences__accession',
        'annotation__component_sequences__description',
        'annotation__component_sequences__organism'
    ]

    acc_id_set = set()
    for row in anno_matches:
        acc_id_set.add(row['annotation__component_sequences__accession'])

    acc_go_edges = getGoIdsByAcc(list(acc_id_set))

    #go_id_set = set()
    #for g in acc_go_edges.values():
    #    go_id_set.update(g)
    #go_nodes = getGoNodes(list(go_id_set))
    go_nodes = getGoNodes(acc_go_edges['ALL'])

    # A "flat" lookup dict, so the GO chart is complete
    # (relying on edges data only gives direct edges)
    acc_go_lookup = dict()
    go_acc_lookup = dict()

    for a in acc_id_set:
        temp_set = set()
        for g in acc_go_edges[a]:
            temp_set.add(g) # add the GO term to the set
            # add any parent nodes (except the root term)
            g_parent = go_nodes[g]["parent_go_id"];
            while g_parent is not None and go_nodes[g_parent]["class_level"] >= 1:
                temp_set.add(g_parent)
                g_parent = go_nodes[g_parent]["parent_go_id"]
        acc_go_lookup[a] = list(temp_set)
        acc_go_lookup[a].sort()

    for a, glist in acc_go_lookup.items():
        for g in glist:
            if g not in go_acc_lookup:
                go_acc_lookup[g] = list()
            go_acc_lookup[g].append(a)

    # More lookup tables to make life easier...
    go_child_lookup = dict()
    for g, go_data in go_nodes.items():
        gp = go_data['parent_go_id']
        if gp not in go_child_lookup:
            go_child_lookup[gp] = list()
        go_child_lookup[gp].append(g)

    tgt_cmp_lookup = dict()
    for row in anno_matches:
        chembl_id = row["annotation__chembl_id_lookup__chembl_id"]
        acc_id = row["annotation__component_sequences__accession"]
        if acc_id not in tgt_cmp_lookup:
            tgt_cmp_lookup[acc_id] = list()
        tgt_cmp_lookup[acc_id].append(chembl_id)

    # Reversed (inversed?) lookup of acc_go_edges
    # TODO: deprecate. replace with go_acc_lookup
    #go_acc_edges = dict()
    #for a, glist in acc_go_edges.items():
    #    if a == 'ALL':
    #        continue
    #    for g in glist:
    #        if g not in go_acc_edges:
    #            go_acc_edges[g] = list()
    #        go_acc_edges[g].append(a)

    # Count number of connections for compounds and targets for initial values
    # of filters.
    ct_max_dict = dict()
    tc_max_dict = dict()
    for m in anno_matches:
        cid = m['annotation__chembl_id_lookup__chembl_id']
        tid = m['annotation__component_sequences__accession']
        if cid in ct_max_dict:
            ct_max_dict[cid] = ct_max_dict[cid] + 1
        else:
            ct_max_dict[cid] = 1
        if tid in tc_max_dict:
            tc_max_dict[tid] = tc_max_dict[tid] + 1
        else:
            tc_max_dict[tid] = 1

    # Get absolute max values
    ct_max = 0
    tc_max = 0
    for v in ct_max_dict.values():
        if v > ct_max:
            ct_max = v
    for v in tc_max_dict.values():
        if v > tc_max:
            tc_max = v

    # Rearrange table info into dict form for easier access
    table_info_dict = dict()
    for i in myAnnoSearch.table_info:
        table_info_dict[i['id']] = i

    context = {
        'ids': ids,
        'anno_info': anno_info,
        'anno_matches': anno_matches,
        'cols': cols,
        'cols_basic': cols_basic,
        'acc_go_edges': acc_go_edges,
        #'go_acc_edges': go_acc_edges,
        'go_nodes': go_nodes,
        'acc_go_lookup': acc_go_lookup,
        'go_acc_lookup': go_acc_lookup,
        'go_child_lookup': go_child_lookup,
        'tgt_cmp_lookup': tgt_cmp_lookup,
        'ct_max_dict': ct_max_dict,
        'tc_max_dict': tc_max_dict,
        'ct_max': ct_max,
        'tc_max': tc_max,
        'table_info_dict': table_info_dict
    }

    return context

def tsFilter2(request):
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    #else:
    #    ids = None

    context = tsFilter2_debug(ids)

    return render(request, 'targetsearch/filter2.html', context)
