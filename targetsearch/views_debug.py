# views_debug.py - Put 'debug' versions of view functions here

#from .helpers import *
from .chembl_search import *
from django.template import Context, Engine

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

    # Gene Ontology stuff starts here

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

    # Gene Ontology stuff ends here
    # Reactome stuff begins here

    acc_reactome_data = getReactomeDataByAcc(list(acc_id_set))
    acc_reactome_edges = getReactomeIdsFromResults(acc_reactome_data)
    reactome_nodes = getReactomeNodesFromResults(acc_reactome_data)

    reactome_acc_lookup = dict()
    for a, rlist in acc_reactome_edges.items():
        for r in rlist:
            if r not in reactome_acc_lookup:
                reactome_acc_lookup[r] = list()
            reactome_acc_lookup[r].append(a)

    # Reactome stuff ends here

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

    org_set = set()
    for a in anno_matches:
        org_set.add(a['annotation__component_sequences__organism'])

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
        'acc_reactome_edges': acc_reactome_edges,
        'reactome_nodes': reactome_nodes,
        'reactome_acc_lookup': reactome_acc_lookup,
        'tgt_cmp_lookup': tgt_cmp_lookup,
        'ct_max_dict': ct_max_dict,
        'tc_max_dict': tc_max_dict,
        'ct_max': ct_max,
        'tc_max': tc_max,
        'table_info_dict': table_info_dict,
        'org_list': list(org_set)
    }

    return context

def tsAnnoFilter1_debug(ids):
    # Get annotation data, just as before
    anno_search = AnnotationWithDrugIndSearch('compound', ids)
    anno_info = anno_search.table_info
    anno_matches = anno_search.get_results()

    if len(anno_matches) == 0:
        context = {
            'ids': ids,
            'anno_info': anno_info,
            'anno_matches': anno_matches,
            'message': 'Search returned no annotation data.'
        }
        return context

    # Rearrange table info into dict form for easier access
    table_info_dict = dict()
    for i in anno_info:
        table_info_dict[i['id']] = i

    # Remove show_extanno. We don't need it here, plus it can't be JSONed
    for m in anno_matches:
        m.pop('annotation__extanno__show_extanno')

    cols = [ k for k in anno_matches[0].keys() ]

    # Minimal set of columns for a compact table
    cols_basic = [
        'annotation__chembl_id_lookup__chembl_id',
        'annotation__molecule_dictionary__pref_name',
        'annotation__drug_mechanism__mechanism_of_action',
        'annotation__component_sequences__accession',
        'annotation__component_sequences__description',
        'annotation__component_sequences__organism'
    ]

    # Combined anno_matches loop
    acc_id_set = set()
    tgt_cmp_lookup = dict()
    ct_max_dict = dict()
    tc_max_dict = dict()
    org_set = set()

    for row in anno_matches:
        acc_id_set.add(row['annotation__component_sequences__accession'])

        chembl_id = row["annotation__chembl_id_lookup__chembl_id"]
        acc_id = row["annotation__component_sequences__accession"]
        if acc_id not in tgt_cmp_lookup:
            tgt_cmp_lookup[acc_id] = list()
        tgt_cmp_lookup[acc_id].append(chembl_id)

        cid = row['annotation__chembl_id_lookup__chembl_id']
        tid = row['annotation__component_sequences__accession']
        if cid in ct_max_dict:
            ct_max_dict[cid] = ct_max_dict[cid] + 1
        else:
            ct_max_dict[cid] = 1
        if tid in tc_max_dict:
            tc_max_dict[tid] = tc_max_dict[tid] + 1
        else:
            tc_max_dict[tid] = 1

        org_set.add(row['annotation__component_sequences__organism'])

    # Get absolute max values
    ct_max = 0
    tc_max = 0
    for v in ct_max_dict.values():
        if v > ct_max:
            ct_max = v
    for v in tc_max_dict.values():
        if v > tc_max:
            tc_max = v

    # Rearrange table info into dict form for easier access TODO: duplicate?
    #table_info_dict = dict()
    #for i in anno_search.table_info:
    #    table_info_dict[i['id']] = i

    # Gene Ontology stuff starts here

    acc_go_edges = getGoIdsByAccFlat(list(acc_id_set), remove_roots=True, add_missing_links=True)
    go_nodes = getGoNodes(acc_go_edges['ALL'])

    go_acc_lookup = dict()
    for a, glist in acc_go_edges.items():
        for g in glist:
            if g not in go_acc_lookup:
                go_acc_lookup[g] = list()
            go_acc_lookup[g].append(a)

    # Gene Ontology stuff ends here

    # Reactome stuff begins here

    acc_reactome_data = getReactomeDataByAcc(list(acc_id_set))
    acc_reactome_edges = getReactomeIdsFromResults(acc_reactome_data)
    reactome_nodes = getReactomeNodesFromResults(acc_reactome_data)

    reactome_acc_lookup = dict()
    for a, rlist in acc_reactome_edges.items():
        for r in rlist:
            if r not in reactome_acc_lookup:
                reactome_acc_lookup[r] = list()
            reactome_acc_lookup[r].append(a)

    # Reactome stuff ends here

    context = {
        'ids': ids,
        'anno_info': anno_info,
        'anno_matches': anno_matches,
        'cols': cols,
        'cols_basic': cols_basic,
        'acc_go_edges': acc_go_edges,
        'go_nodes': go_nodes,
        'go_acc_lookup': go_acc_lookup,
        'acc_reactome_edges': acc_reactome_edges,
        'reactome_nodes': reactome_nodes,
        'reactome_acc_lookup': reactome_acc_lookup,
        'tgt_cmp_lookup': tgt_cmp_lookup,
        #'ct_max_dict': ct_max_dict,
        #'tc_max_dict': tc_max_dict,
        'ct_max': ct_max,
        'tc_max': tc_max,
        'table_info_dict': table_info_dict,
        'org_list': list(org_set)
    }

    return context

def tsAnnoFilter2_debug(ids):
    go_acc_lookup = getAccIdsByGo1(ids)

    # Combined go_acc_lookup loop
    acc_id_set = set()
    acc_go_edges = dict()
    for go_id, acc_list in go_acc_lookup.items():
        acc_id_set.update(acc_list)

        for a in acc_list:
            if a not in acc_go_edges:
                acc_go_edges[a] = [go_id]
            else:
                acc_go_edges[a].append(go_id)
    acc_go_edges['ALL'] = list(go_acc_lookup.keys())

    anno_search = AnnotationWithDrugIndSearch('target', list(acc_id_set))
    anno_info = anno_search.table_info
    anno_matches = anno_search.get_results()

    if len(go_acc_lookup) == 0:
        context = {
            'ids': ids,
            'anno_info': anno_info,
            'anno_matches': anno_matches,
            'message': 'GO search returned no Accession ID data.'
        }
        return context
    elif len(anno_matches) == 0:
        context = {
            'ids': ids,
            'anno_info': anno_info,
            'anno_matches': anno_matches,
            'message': 'Accession ID search returned no compound data.'
        }
        return context

    table_info_dict = dict()
    for i in anno_info:
        table_info_dict[i['id']] = i

    for m in anno_matches:
        m.pop('annotation__extanno__show_extanno')

    cols = [ k for k in anno_matches[0].keys() ]

    cols_basic = [
        'annotation__component_sequences__accession',
        'annotation__component_sequences__description',
        'annotation__component_sequences__organism',
        'annotation__chembl_id_lookup__chembl_id',
        'annotation__molecule_dictionary__pref_name',
        'annotation__drug_mechanism__mechanism_of_action'
    ]

    # Combined anno_matches loop
    tgt_cmp_lookup = dict()
    ct_max_dict = dict()
    tc_max_dict = dict()
    org_set = set()

    for row in anno_matches:
        chembl_id = row["annotation__chembl_id_lookup__chembl_id"]
        acc_id = row["annotation__component_sequences__accession"]
        if acc_id not in tgt_cmp_lookup:
            tgt_cmp_lookup[acc_id] = list()
        tgt_cmp_lookup[acc_id].append(chembl_id)

        cid = row['annotation__chembl_id_lookup__chembl_id']
        tid = row['annotation__component_sequences__accession']
        if cid in ct_max_dict:
            ct_max_dict[cid] = ct_max_dict[cid] + 1
        else:
            ct_max_dict[cid] = 1
        if tid in tc_max_dict:
            tc_max_dict[tid] = tc_max_dict[tid] + 1
        else:
            tc_max_dict[tid] = 1

        org_set.add(row['annotation__component_sequences__organism'])

    # Get absolute max values
    ct_max = 0
    tc_max = 0
    for v in ct_max_dict.values():
        if v > ct_max:
            ct_max = v
    for v in tc_max_dict.values():
        if v > tc_max:
            tc_max = v

    # Gene Ontology stuff starts here

    # acc_go_edges = see combined loop above
    go_nodes = getGoNodes(acc_go_edges['ALL'])

    # Gene Ontology stuff ends here

    # Reactome stuff starts here

    # The template expects this variable, so fill it with empty lists for now
    acc_reactome_edges = dict()
    for acc in acc_id_set:
        acc_reactome_edges[acc] = list()
    acc_reactome_edges['ALL'] = list()

    # Reactome stuff ends here

    context = {
        'ids': ids,
        'anno_info': anno_info,
        'anno_matches': anno_matches,
        'cols': cols,
        'cols_basic': cols_basic,
        'acc_go_edges': acc_go_edges,
        'go_nodes': go_nodes,
        'go_acc_lookup': go_acc_lookup,
        'acc_reactome_edges': acc_reactome_edges,
        'reactome_nodes': dict(),
        'reactome_acc_lookup': dict(),
        'tgt_cmp_lookup': tgt_cmp_lookup,
        'ct_max': ct_max,
        'tc_max': tc_max,
        'table_info_dict': table_info_dict,
        'org_list': list(org_set)
    }

    return context
