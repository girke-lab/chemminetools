# helpers.py - Anything that doesn't quite fit in views or views_debug
# probably should go here

import psycopg2
import psycopg2.sql as sql
import requests
import json
import os
import pybel
import tempfile
import subprocess
import shutil
import pprint
import sqlite3
import csv
from psycopg2.extras import NamedTupleCursor
from django.conf import settings
from collections import OrderedDict
#from .chembl_search import *

def ts_paralog_cache():
    try:
        return settings.TS_PARALOG_CACHE
    except AttributeError as e:
        return False

# Useful to avoid waiting for Ensembl to respond during development
if ts_paralog_cache():
    from sqlitedict import SqliteDict

def groupBy(keyFn, row_data):
    """Group rows in row_data by the result of keyFn.

    Arguments:
    keyFn -- A function that returns the cell data of the column to group by.
             Needs to be tailored to the structure of row_data.
    row_data -- A table organized as a list of row data structures.

    Returns:
    { 'key_id1': [ row1, row2, row3, ...],
      'key_id2': [ row4, row5, ... ],
      ...}
    """
    temp_dict = {}
    for row in row_data:
        key = keyFn(row) # Get key/cell data corresponding to the groupBy column for this row
        a = temp_dict.get(key,[]) # Get the sub-table for whatever we have so far for the given key
        a.append(row) # Attach this row to the grouped sub-table
        temp_dict[key] = a # Put sub-table back into temp_dict

    return temp_dict

def runQuery(query, values, cursor_factory=psycopg2.extras.NamedTupleCursor):
    """Queries the local ChEMBL database.

    Currently assumes a PostgreSQL database, so write queries according to psycopg2
    formats. A lot of old code assumes results are returned with NamedTupleCursor
    enabled, hence it's the default for CMT. New code should adopt DictCursor.

    Arguments:
    query -- PostgreSQL query. Please be cognizant of SQL injection attacks.
    values -- Values to fill into the query.
    cursor_factory -- psycopg2 cursor factory. Pass the cursor object directly
        (psycopg2.extras.X). Alternatively, set to None for default psycopg2 behavior.

    Returns a list of records. The record structure will depend on the cursor factory.
    """
    dbhost = settings.CHEMBL_DB['DBHOST']
    dbname = settings.CHEMBL_DB['DBNAME']
    dbuser = settings.CHEMBL_DB['DBUSER']
    dbpass = settings.CHEMBL_DB['DBPASS']

    with psycopg2.connect(host=dbhost, dbname=dbname, user=dbuser, password=dbpass) as conn:
        if cursor_factory is None:
            cur = conn.cursor()
        else:
            cur = conn.cursor(cursor_factory=cursor_factory)
        cur.execute(query, values)
        return cur.fetchall()

def get_chembl_sdfs(chemblIds):
    data = runQuery("""SELECT chembl_id||molfile, pref_name
      FROM molecule_dictionary JOIN
           compound_structures USING(molregno)
      WHERE molecule_dictionary.chembl_id IN %s""",(tuple(chemblIds),))

    sdf_list = list()
    # Attach pref_name, if it exists
    for row in data:
        molfile = row[0]
        if row.pref_name is not None:
            molfile += ('\n> <NAME>\n' + row.pref_name)
        sdf_list.append(molfile)

    return sdf_list

def get_chembl_smiles(chemblIds):
    data = runQuery("""
      SELECT canonical_smiles||' '||chembl_id
      FROM chembl_id_lookup JOIN
           compound_structures ON(entity_id=molregno)
      WHERE chembl_id_lookup.chembl_id IN %s""",(tuple(chemblIds),))
    return [row[0] for row in data]

def get_chembl_ids_with_struct(ids):
    data = runQuery("""
        SELECT chembl_id
        FROM molecule_dictionary
        JOIN compound_structures USING(molregno)
        WHERE molecule_dictionary.chembl_id IN %s""", (tuple(ids),))
    return [ d.chembl_id for d in data ]

def batchQuery(runQuery,ids, batchSize=1000):
    finalResult = []
    for i in range(0,len(ids),batchSize):
        batchResult = runQuery(ids[i:i+batchSize])
        finalResult = finalResult + batchResult
    return finalResult

def mapToChembl(unknownIds, sourceId):
    chemblIds = {}
    for unknownId in unknownIds:
        req = requests.get("https://www.ebi.ac.uk/unichem/rest/src_compound_id/"+unknownId+"/"+str(sourceId)+"/1")
        #print("\n\nreq.json: "+str(req.json()))
        result = req.json()
        #print("got result: "+str(result))
        if "error" in result:
            raise Exception(result["error"])
        if isinstance(result,list) and len(result) > 0 and ("src_compound_id" in result[0]) :
            chemblIds[result[0]["src_compound_id"] ] = unknownId

    #print("final chembl ids: "+str(chemblIds))

    return chemblIds

def mapToUniprot(unknownIds, sourceId):
    return uniprotMap(unknownIds, sourceId, "ACC")

def uniprotMap(unknownIds, a, b, reverse=True):
    """Generic UniProt ID mapping function. Translates a list of unknownIds
    from type a to type b."""
    url = "https://www.uniprot.org/uploadlists/"
    params = {
        "from": a,
        "to": b,
        "format": "tab",
        "query": " ".join(unknownIds)
    }

    results = requests.get(url, params=params).text
    results_dict = {}
    try:
        for line in results.splitlines():
            if "From\tTo" not in line:
                row = line.split("\t")
                if reverse:
                    results_dict[row[1]] = row[0]
                else:
                    results_dict[row[0]] = row[1]
    except Exception as e:
        print("Oops, something happened. Take a look:")
        print(results)
        raise e

    return results_dict

def getUniChemSources():
    sources = {}
    sourceReq = requests.get("https://www.ebi.ac.uk/unichem/rest/src_ids/")
    sourceIds = [ source["src_id"] for source in sourceReq.json()]

    for sourceId in sourceIds:
        infoReq = requests.get(" https://www.ebi.ac.uk/unichem/rest/sources/"+sourceId)
        data = infoReq.json()[0]
        #print(str(data["src_id"])+"\t"+data["name_label"])
        sources[data["src_id"]] = data["name_label"]

    return sources;

def compoundNameAutocomplete(nameQuery):
    sqlQuery = sql.SQL("""SELECT DISTINCT chembl_id, synonyms
                          FROM molecule_synonyms
                          JOIN molecule_dictionary USING(molregno)
                          WHERE synonyms ILIKE %s
                          ORDER BY synonyms""")

    return runQuery(sqlQuery, ("%"+nameQuery.replace(' ','%')+"%",))

def targetNameAutocomplete(nameQuery):
    sqlQuery = sql.SQL("""SELECT accession, description, organism
                          FROM component_sequences
                          WHERE description||' '||organism ILIKE %s
                          ORDER BY description, organism""")

    return runQuery(sqlQuery, ('%'+nameQuery.replace(' ','%')+'%',))
    
def compoundSummaryInfo(chemblId):
    result = runQuery("""SELECT full_mwt,full_molformula, pref_name,chembl_id,molecule_type 
                    FROM  molecule_dictionary JOIN compound_properties USING(molregno)  
                    WHERE chembl_id=%s""",(chemblId,))

    if len(result) == 1:
        return result[0];
    return None

def targetSummaryInfo(accession):
    sqlQuery = sql.SQL(""" select accession, description, organism, 
                            (select string_agg(component_synonym,',') 
                                from component_synonyms as cs 
                                where cs.component_id=cseq.component_id 
                                      and syn_type='GENE_SYMBOL'
                            ) as gene_names 
                           from component_sequences as cseq 
                           where accession = %s """)

    result = runQuery(sqlQuery, (accession,))
    if len(result) == 1:
        return result[0];
    return None


def getEntityType(entityId):
    """when we have an ID, but don't know whether its a CHEMBL ID
       or an accession or something else
    """
    chemblTest = sql.SQL("""SELECT DISTINCT 1 
                            FROM molecule_dictionary
                            WHERE chembl_id = %s""")
    accessionTest = sql.SQL("""SELECT DISTINCT 1 
                               FROM component_sequences
                               WHERE accession = %s""")
    if len(runQuery(chemblTest,(entityId,)))== 1:
        return "compound"
    elif len(runQuery(accessionTest,(entityId,))) == 1:
        return "target"
    else:
        return "unknown"

def getChemblPNG(chembl_id, size=500, mwt_limit=None, write_file=False,
                 filename=None, shrink=True):
    """Given a ChEMBL ID, generate a PNG of the molecular structure with
    OpenBabel. Return either the PNG bytes, or a string of the filename.

    Arguments:
    chembl_id (str) -- ChEMBL ID of compound
    size (int, str) -- Dimension of PNG file in pixels ('p' option in OB).
                       Default 500
    mwt_limit (int, float, None) -- If the molecular weight of the compound is
                             greater than this limit, raise an exception.
                             Default None.
    write_file (bool) -- If True, write the PNG to a file and return the file
                         path instead of the PNG bytes. An existing file will
                         be overwritten. Default False.
    filename (str) -- Path to output file. If None, create a temporary file.
                      Default None.
    shrink (bool) -- If true, run a utility (currently OptiPNG) to shrink the
                     generated PNG. Default True.
    """
    data = runQuery("""SELECT standard_inchi, full_mwt
                       FROM molecule_dictionary
                       JOIN compound_properties USING(molregno)
                       JOIN compound_structures USING(molregno)
                       WHERE chembl_id = %s
                       LIMIT 1""", (chembl_id,))
    inchi = data[0].standard_inchi
    mwt = float(data[0].full_mwt)

    if mwt_limit is not None and mwt > mwt_limit:
        raise Exception("{} exceeded mwt_limit (mwt={}, mwt_limit={})".format(chembl_id, mwt, mwt_limit))

    mymol = pybel.readstring('inchi', inchi)

    fp = tempfile.NamedTemporaryFile(delete=False)
    fp.close()
    tmpfile = fp.name

    mymol.write(format='_png2', filename=tmpfile, overwrite=True, opt={'p': size})

    if shrink:
        cmd = ['optipng', tmpfile] # Losslessly shrinks PNG by ~15-25%
        subprocess.run(cmd, check=True)

    if write_file:
        if filename is not None:
            return shutil.move(tmpfile, filename)
        else:
            return tmpfile
    else:
        with open(tmpfile, 'rb') as f:
            img = f.read()
        os.remove(tmpfile)
        return img

def getChemblSVG(chembl_id, mwt_limit=None, filename=None):
    """Given a ChEMBL ID, generate a PNG of the molecular structure with
    OpenBabel. Return either the SVG data string, or a string of the filename.

    Arguments:
    chembl_id (str) -- ChEMBL ID of compound
    mwt_limit (int, float, None) -- If the molecular weight of the compound is
                             greater than this limit, raise an exception.
                             Default None.
    filename (str) -- Path to output file. If None, return SVG data.
                      Default None."""

    data = runQuery("""SELECT molfile, full_mwt
                       FROM molecule_dictionary
                       JOIN compound_properties USING(molregno)
                       JOIN compound_structures USING(molregno)
                       WHERE chembl_id = %s
                       LIMIT 1""", (chembl_id,))
    if len(data) == 0:
        raise Exception("ChEMBL ID not found for: {}".format(chembl_id))
    mdl = data[0].molfile
    mwt = float(data[0].full_mwt)

    if mwt_limit is not None and mwt > mwt_limit:
        raise Exception("{} exceeded mwt_limit (mwt={}, mwt_limit={})".format(chembl_id, mwt, mwt_limit))

    mymol = pybel.readstring('mdl', mdl)
    svg = mymol.write(format='svg', opt={'d': ''})

    if filename is None:
        return svg
    else:
        with open(filename, 'w') as f:
            f.write(svg)
        return filename

def getTargetDesc(accList):
    data = runQuery("""SELECT accession, description, organism
                          FROM component_sequences
                          WHERE accession IN %s""", (tuple(accList),))

    descDict = dict()
    for r in data:
        descDict[r.accession] = "{} ({})".format(r.description, r.organism)

    return descDict

def homologsAjax(ids, homolog_type):
    if len(ids) == 0:
        raise Exception('Empty list of "ids".')

    if homolog_type == 'paralog':
        r = getParalogs(ids)
        paralogs = r['paralogs']
        uniprotIdMap = r['uniprotIdMap']
        idErrors = set(ids) - set(paralogs.keys())
        cacheMisses = None
        uids = gatherUniprotIds(paralogs)
        descriptions = getTargetDesc(uids)
        chemblErrors = set(uids) - set(descriptions.keys())
    elif homolog_type == 'paralog-cache':
        r = getParalogsCached(ids)
        paralogs = r['paralogs']
        cacheMisses = r['cacheMisses']
        uniprotIdMap = r['uniprotIdMap']
        idErrors = set(ids) - set(paralogs.keys())
        uids = gatherUniprotIds(paralogs)
        descriptions = getTargetDesc(uids)
        chemblErrors = set(uids) - set(descriptions.keys())
    else:
        raise Exception('Invalid value for "homolog_type".')

    ajaxResponse = {'success': True,
                    'paralogs': paralogs,
                    'idErrors': list(idErrors),
                    'cacheMisses': cacheMisses,
                    'uniprotIdMap': uniprotIdMap,
                    'descriptions': descriptions,
                    'chemblErrors': list(chemblErrors),
                    }
    return ajaxResponse

def gatherUniprotIds(paralogsDict):
    ids = set()
    for k, v in paralogsDict.items():
        ids.add(k)
        for p in v:
            ids.add(p['uniprot'])
    return ids

def getParalogs(queryIds):
    ids = set(queryIds)
    base = "http://www.ensembl.org/biomart/martservice?query="
    uniprotQuery= """<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
            <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                <Filter name = "uniprot_gn_id" value ="{}" />
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "hsapiens_paralog_ensembl_gene" />
                <Attribute name = "hsapiens_paralog_perc_id" />
                <Attribute name = "hsapiens_paralog_perc_id_r1" />
                <Attribute name = "hsapiens_paralog_associated_gene_name" />
                <Attribute name = "external_gene_name" />
           </Dataset>
        </Query>""".format(",".join(ids)).replace("\n","")
                              #<Attribute name = "hsapiens_paralog_ds" />
                #<Attribute name = "hsapiens_paralog_dn" />

    paralogResults = requests.get(base + uniprotQuery).text

    # Stop early if there are no results
    if paralogResults == '':
        return {}
    if 'ERROR' in paralogResults:
        raise Exception(paralogResults)

    paralogs = {}
    geneIds = set()
#    print("paralog results: \n"+str(paralogResults))

    for line in paralogResults.splitlines():
        row = line.split("\t")
        geneIds.add(row[0])
        geneIds.add(row[1])
        if row[0] not in paralogs:
            paralogs[row[0]] = list()
        d = dict()
        d["ensembl_gene"] = row[1]
        d["perc_id"] = row[2]
        d["perc_id_r1"] = row[3]
        d["gene_name"] = row[4]
        paralogs[row[0]].append(d)

        #paralogs[row[1]]=row[0]

#    print("gene ids: "+str(geneIds))
#    print("paralogs: \n"+str(paralogs))

    ensemblQuery= """<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
            <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                <Filter name = "ensembl_gene_id" value ="{}"/>
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "uniprot_gn_symbol" />
                <Attribute name = "uniprotswissprot" />
                <Attribute name = "uniprotsptrembl" />
            </Dataset>
        </Query>""".format(",".join(geneIds)).replace("\n","")


    ensemblResults = requests.get(base + ensemblQuery).text

    if 'ERROR' in ensemblResults:
        raise Exception(ensemblResults)

    geneToUniprot = {}
    uniprotIdMap = {}
    for line in ensemblResults.splitlines():
        row = line.split("\t")
        #print(row)
        ensembl_id = row[0]
        gene_name = row[1]
        uniprot_id = row[2]
        trembl_id = row[3]
        if ensembl_id not in geneToUniprot:
            geneToUniprot[ensembl_id] = {"uniprot":"","trembl":[]}
        if uniprot_id != "":
            geneToUniprot[ensembl_id]["uniprot"] = uniprot_id
            uniprotIdMap[uniprot_id] = {"ensembl": ensembl_id,
                                        "gene_name": gene_name}
        elif trembl_id != "":
            geneToUniprot[ensembl_id]["trembl"] += [trembl_id]

    #print("geneToUniprot: "+str(geneToUniprot))

    results = {}
    for geneId in paralogs:
        hits = paralogs[geneId]
        id = geneId
        if geneId in geneToUniprot: #translate id if we can
            id = geneToUniprot[id]["uniprot"]

        for hit in hits:  # for each hit, translate the hit gene name if we can
            e = hit["ensembl_gene"]
            if e in geneToUniprot:
                hit["trembl"] = geneToUniprot[e]["trembl"]
                hit["uniprot"] = geneToUniprot[e]["uniprot"]

        results[id] = hits


    #print("results: \n")
    #pp = pprint.PrettyPrinter(indent=4)
    #pp.pprint(results)
    return {"paralogs": results, "uniprotIdMap": uniprotIdMap}

#getParalogs(['P35354','P23219','P02768'])

def getParalogsCached(queryIds):
    with SqliteDict('/srv/getParalogs/paralogs.sqlite') as paralogsCache,\
         SqliteDict('/srv/getParalogs/idmap.sqlite') as idMapCache:
        querySet = set(queryIds) - set(paralogsCache.keys())
        if len(querySet) != 0:
            #print("Cache misses: {}".format(" ".join(querySet)))
            gp_results = getParalogs(querySet)
            gp_paralogs = gp_results['paralogs']
            gp_idmap = gp_results['uniprotIdMap']

            for k, v in gp_paralogs.items():
                paralogsCache[k] = v

            for k, v in gp_idmap.items():
                idMapCache[k] = v

            # Cache query errors as None
            missingIds = querySet - set(gp_paralogs.keys())
            for i in missingIds:
                paralogsCache[i] = None

            paralogsCache.commit()
            idMapCache.commit()

        paralogs = dict()
        for q in queryIds:
            r = paralogsCache.get(q)
            if r is not None:
                paralogs[q] = r

        uniprotIdMap = dict()
        #for k, v in paralogs.items():
        #    uniprotIdMap[k] = idMapCache[k]
            #for p in v:
            #    uniprotIdMap[p['uniprot']] = idMapCache[p['uniprot']]
        for i in gatherUniprotIds(paralogs):
            uniprotIdMap[i] = idMapCache[i]

        return {'paralogs': paralogs,
                'cacheMisses': list(querySet),
                'uniprotIdMap': uniprotIdMap}

def queryExtAnno(chembl_id):
    """Retrieves all available extended annotations for a given ChEMBL ID.

    Currently not used anywhere, and is slow on SQLite. It's better to use
    checkExtAnno() and getExtAnno(). Consider deprecating."""
    with sqlite3.connect(settings.ANNO_DB) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        query = """SELECT chembl_id, drugbank_id, drugage_id, cmap_id, lincs_id
            FROM id_mapping
            WHERE chembl_id = ?"""
        cur.execute(query, (chembl_id,))
        id_map = cur.fetchone()

        if id_map is None:
            return None

        results = dict()
        for a in [("drugbank_id","DrugBankAnnot"),
                  ("drugage_id","drugAgeAnnot"),
                  ("cmap_id","cmapAnnot"),
                  ("lincs_id","lincsAnnot")]:
            key, table = a
            if id_map[key] is None:
                results[key] = None
            else:
                query = """SELECT * FROM {table} WHERE {col} = ?""".format(table=table, col=key)
                d = dict()
                cur.execute(query, (id_map[key],))
                anno_results = cur.fetchone()
                for k in anno_results.keys():
                    d[k] = anno_results[k]
                results[key] = d

        return results

def checkExtAnno(chembl_id):
    """Checks the availability of an extended annotation for a given ChEMBL ID.

    Returns a dict object with translated IDs, or None if a translation DNE."""
    with sqlite3.connect(settings.ANNO_DB) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        query = """SELECT drugbank_id as drugbank, drugage_id as drugage, cmap_id as cmap, lincs_id as lincs
                   FROM id_mapping
                   WHERE chembl_id = ?"""
        cur.execute(query, (chembl_id,))
        id_map = cur.fetchone()

        if id_map is None:
            return None

        results = dict()
        results.update(id_map)

        return results

def getExtAnno(chembl_id, db):
    with sqlite3.connect(settings.ANNO_DB) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        if db == "drugbank":
            join_col = "drugbank_id"
            join_table = "DrugBankAnnot"
        elif db == "drugage":
            join_col = "drugage_id"
            join_table = "drugAgeAnnot"
        elif db == "cmap":
            join_col = "cmap_id"
            join_table = "cmapAnnot"
        elif db == "lincs":
            join_col = "lincs_id"
            join_table = "lincsAnnot"
        else:
            raise Exception("Invalid db value: {}".format(db))
        query = """SELECT id_mapping.chembl_id, {join_table}.*
                   FROM id_mapping
                   JOIN {join_table} using({join_col})
                   WHERE chembl_id = ?""".format(join_col=join_col, join_table=join_table)
        cur.execute(query, (chembl_id,))
        anno_results = cur.fetchone()

    if anno_results is None:
        return None

    results = dict()
    results.update(anno_results)

    return results

def getExtAnno2(id_type, cid):
    if id_type == "drugbank_id":
        table = "DrugBankAnnot"
    elif id_type == "drugage_id":
        table = "drugAgeAnnot"
    elif id_type == "cmap_id":
        table = "cmapAnnot"
    elif id_type == "lincs_id":
        table = "lincsAnnot"
    else:
        raise Exception("Invalid id_type value: {}".format(id_type))

    with sqlite3.connect(settings.ANNO_DB) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        query = """SELECT *
                   FROM {table}
                   WHERE {id_type} = ?""".format(table=table, id_type=id_type)
        cur.execute(query, (cid,))
        anno_results = cur.fetchone()

    if anno_results is None:
        return None
    else:
        results = dict()
        results.update(anno_results)
        return results

def getGoData(acc_id_list):
    """Given a list of UniProt Accession IDs, query the ChEMBL database for
    Gene Ontology classifications.

    ['P23219', 'P35354'] => {   'P23219' : [r1,r2,...],
                                'P35354' : [r1,r2,...]}
    """
    go_data_dict = dict()
    for acc in acc_id_list:
        data = runQuery(""" SELECT  component_sequences.accession,
                                    component_sequences.component_id,
                                    component_go.go_id,
                                    go_classification.parent_go_id,
                                    go_classification.pref_name,
                                    go_classification.class_level,
                                    go_classification.aspect,
                                    go_classification.path
                            FROM component_sequences
                            JOIN component_go USING(component_id)
                            JOIN go_classification USING(go_id)
                            WHERE accession = %s
                            """, (acc,), cursor_factory=psycopg2.extras.DictCursor)
        go_data_dict[acc] = [dict(row) for row in data]
    return go_data_dict

def getGoIdsByAcc_OLD(acc_id_list):
    """Given a list of UniProt Accession IDs, query the ChEMBL database for
    associated Gene Ontology IDs.

    For the Accession IDs, only the most specific GO IDs (leaf nodes) will be
    connected. This is to make it easier to build a proper GO tree structure.

    In addition, there is a special list under the "ALL" key, which contains
    all GO IDs from the original query. This list should be fed into getGoNodes
    to obtain further data.

    Example:

    ['P23219', 'P35354'] => {   'P23219' : [go_id,go_id,...],
                                'P35354' : [go_id,go_id,...]}"""
    acc_go_dict = dict()
    all_go_ids = set()
    for acc in acc_id_list:
        data = runQuery(""" SELECT  component_sequences.accession,
                                    component_go.go_id,
                                    go_classification.parent_go_id
                            FROM component_sequences
                            JOIN component_go USING(component_id)
                            JOIN go_classification USING(go_id)
                            WHERE accession = %s
                            """, (acc,), cursor_factory=psycopg2.extras.DictCursor)
        # Determine which nodes are branches (internal) and which are leaves (external).
        # Assuming it might be possible for a GO node to be both a branch or a leaf depending
        # on specific drug-target context, we make this determination for each drug-target,
        # rather than at the greater GO tree context.
        bnodes = [row['parent_go_id'] for row in data if row['parent_go_id'] is not None]
        #bnodes = [row['parent_go_id'] for row in data]
        #acc_go_dict[acc] = [row['go_id'] for row in data if row['go_id'] not in bnodes]
        go_id_list = list()
        for row in data:
            g = row['go_id']
            all_go_ids.add(g)
            # Do not trust the ChEMBL GO table to provide all necessary nodes for a
            # proper graph. For example, ChEMBL 27 is missing a record linking
            # P23219 to GO:0003824 (parent of GO:0016491, which is linked)
            if row['parent_go_id'] is not None:
                all_go_ids.add(row['parent_go_id'])
            if g not in bnodes:
                go_id_list.append(g)
        acc_go_dict[acc] = go_id_list
    acc_go_dict['ALL'] = list(all_go_ids)
    return acc_go_dict

def getGoIdsByAcc(acc_id_list, flat = False):
    """Given a list of UniProt Accession IDs, query the ChEMBL database for
    associated Gene Ontology IDs.

    If "flat" is False (default), only the most specific GO IDs (leaf nodes)
    will be connected to the Accession ID, so a proper GO tree structure can be
    built later. If True, all GO IDs will be connected directly to the
    Accession ID. NOTE: use getGoIdsByAccFlat() instead.

    In addition, there is a special list under the "ALL" key, which contains
    all GO IDs from the original query. This list should be fed into getGoNodes
    to obtain further data.

    Example:

    ['P23219', 'P35354'] => {   'P23219' : [go_id,go_id,...],
                                'P35354' : [go_id,go_id,...],
                                'ALL' : [go_id,go_id,go_id,...]}"""
    acc_go_dict = dict()
    all_go_ids = set()
    parent_dict = dict()

    # Query the entire go_id -> parent_go_id relation table. It's not that
    # large as of ChEMBL 27
    data1 = runQuery("""SELECT go_id, parent_go_id FROM go_classification""",
                     None, cursor_factory=psycopg2.extras.DictCursor)
    for row in data1:
        parent_dict[row['go_id']] = row['parent_go_id']

    for acc in acc_id_list:
        data2 = runQuery("""SELECT  component_sequences.accession,
                                    component_go.go_id,
                                    go_classification.parent_go_id
                            FROM component_sequences
                            JOIN component_go USING(component_id)
                            JOIN go_classification USING(go_id)
                            WHERE accession = %s""",
                            (acc,), cursor_factory=psycopg2.extras.DictCursor)
        # Determine which nodes are branches (internal) and which are leaves (external).
        # Assuming it might be possible for a GO node to be both a branch or a leaf depending
        # on specific drug-target context, we make this determination for each drug-target,
        # rather than at the greater GO tree context.
        if not flat:
            bnodes = [row['parent_go_id'] for row in data2 if row['parent_go_id'] is not None]

        leaf_list = list()
        for row in data2:
            gid = row['go_id']

            if flat or gid not in bnodes:
                leaf_list.append(gid)
        acc_go_dict[acc] = leaf_list.copy()

        # Manually add parent nodes up to the root, since the component_go
        # table in ChEMBL does not provide all the necessary links for a proper
        # graph. For example, ChEMBL 27 is missing a record linking
        # P23219 to GO:0003824 (parent of GO:0016491, which is linked)
        for l in leaf_list:
            all_go_ids.add(l)
            p = parent_dict[l]
            while p is not None:
                all_go_ids.add(p)
                if flat:
                    acc_go_dict[acc].append(p)
                p = parent_dict[p]
    acc_go_dict['ALL'] = list(all_go_ids)
    return acc_go_dict

def getGoIdsByAccFlat(acc_id_list, remove_roots = True, add_missing_links = True):
    """Given a list of UniProt Accession IDs, query the ChEMBL database for
    associated Gene Ontology IDs.

    Unlike getGoIdsByAcc(), this function connects all GO IDs directly to the
    Accession ID (i.e. no tree structure).

    If "remove_roots" is True (default), the root nodes for the three GO aspects
    will not be included in the result set.

    If "add_missing_links" is True (default), intermediate nodes will be added
    to the result set by climbing up the parent nodes. Otherwise, the GO terms
    from the ChEMBL database will be returned as is.

    In addition, there is a special list under the "ALL" key, which contains
    all GO IDs from the original query. This list should be fed into getGoNodes
    to obtain further data.

    Example:

    ['P23219', 'P35354'] => {   'P23219' : [go_id,go_id,...],
                                'P35354' : [go_id,go_id,...],
                                'ALL' : [go_id,go_id,go_id,...]}"""
    acc_go_dict = dict() # returned result set
    all_go_ids = set()
    parent_dict = dict() # d[go_id] => parent_go_id

    # Query the entire go_id => parent_go_id relation table. It's not that
    # large as of ChEMBL 27
    data1 = runQuery("""SELECT go_id, parent_go_id FROM go_classification""",
                     None, cursor_factory=psycopg2.extras.DictCursor)
    # Build parent_dict lookup. For finding missing links
    for row in data1:
        parent_dict[row['go_id']] = row['parent_go_id']

    for acc in set(acc_id_list):
        #print("processing {}".format(acc))
        data2 = runQuery("""SELECT  component_sequences.accession,
                                    component_go.go_id,
                                    go_classification.parent_go_id
                            FROM component_sequences
                            JOIN component_go USING(component_id)
                            JOIN go_classification USING(go_id)
                            WHERE accession = %s""",
                            (acc,), cursor_factory=psycopg2.extras.DictCursor)

        acc_go_dict[acc] = list()
        for row in data2:
            # Skip root nodes if requested
            if row['parent_go_id'] is None and remove_roots is True:
                pass
            else:
                acc_go_dict[acc].append(row['go_id'])

        # Manually add parent nodes if requested, since the component_go
        # table in ChEMBL does not provide all intermediate links
        # For example, ChEMBL 27 is missing a record linking
        # P23219 to GO:0003824 (parent of GO:0016491, which is linked)
        if add_missing_links:
            for g in acc_go_dict[acc]:
                p = parent_dict[g]
                if parent_dict[p] is None and remove_roots is True:
                    #print("excluding {} -> {}".format(p, parent_dict[p]))
                    pass
                elif p in acc_go_dict[acc]:
                    #print("skipping {}".format(p))
                    pass
                else:
                    #print("appending {}".format(p))
                    acc_go_dict[acc].append(p)

    # Gather all_go_ids
    for v in acc_go_dict.values():
        all_go_ids.update(v)
    acc_go_dict['ALL'] = list(all_go_ids)

    return acc_go_dict

def getAccIdsByGo1(go_id_list):
    """Given a list of Gene Ontology IDs, query the ChEMBL database for
    associated UniProt Accession IDs. Basically the reverse of
    getGoIdsByAccFlat().

    Example:

    ['GO:0005515', 'GO:0016020'] => {   'GO:0005515' : [acc_id,acc_id,...],
                                        'GO:0016020' : [acc_id,acc_id,...]}
    """

    go_acc_dict = dict() # returned result set

    for go_id in set(go_id_list):
        data = runQuery(""" SELECT  component_go.go_id,
                                    component_sequences.accession
                            FROM component_go
                            JOIN component_sequences USING(component_id)
                            WHERE go_id = %s""",
                            (go_id,), cursor_factory=psycopg2.extras.DictCursor)

        go_acc_dict[go_id] = [ row['accession'] for row in data ]

    return go_acc_dict

def getGoNodes(go_id_list):
    """Given a list of Gene Ontology IDs, query the ChEMBL database for
    GO node data. The expected use of this function is to give it the
    ['ALL'] subset of the output of getGoIdsByAcc, but a manually created
    list will be accepted just as well. Note that parent GO nodes are not
    automatically determined by this function (but getGoIdsByAcc does).

    Example:
    GO:0006629 => lipid metabolic process
    GO:0008152 => metabolic process
    GO:0008150 => biological process

    ['GO:0006629', 'GO:0008152', 'GO:0008150'] => { 'GO:0006629' : {data_dict},
                                                    'GO:0008152' : {data_dict},
                                                    'GO:0008150' : {data_dict}}
    """
    go_node_dict = dict()
    data = runQuery(""" SELECT go_id, parent_go_id, pref_name, class_level, aspect, path
                        FROM go_classification
                        WHERE go_id IN %s
                        """, (tuple(go_id_list),), cursor_factory=psycopg2.extras.DictCursor)
    for row in data:
        go_node_dict[row['go_id']] = dict(row)
    return go_node_dict

def getReactomeDataByAcc(acc_id_list):
    """Given a list of UniProt Accession IDs, query the Reactome database for
    associated Reactome pathway IDs and info.

    ['P23219', 'P35354'] => {   'P23219' : [{data_dict}, {data_dict}, ...],
                                'P35354' : [{data_dict}, {data_dict}, ...]}
    """
    results = dict()

    conn = sqlite3.connect(settings.REACTOME_DB)
    conn.row_factory = sqlite3.Row
    query = """SELECT * FROM uniprot WHERE uniprot_id = ?"""
    for acc_id in acc_id_list:
        l = list()
        for row in conn.execute(query, (acc_id,)):
            l.append(dict(row))
        results[acc_id] = l

    return results

def getReactomeNodesFromResults(reactome_results):
    """Given the output of getReactomeDataByAcc, extract data meant for the
    graph display.

    { 'P23219' : [dicts], 'P35354' : [dicts]} => {  'R-HSA-140180' : {data_dict},
                                                    'R-HSA-2162123' : {data_dict},
                                                    ...}
    """
    reactome_node_dict = dict()

    for k1, v1 in reactome_results.items():
        for row in v1:
            if row['reactome_id'] not in reactome_node_dict:
                r = dict()
                for i in ['reactome_id', 'url', 'event_name', 'evidence_code', 'species']:
                    r[i] = row[i]
                reactome_node_dict[row['reactome_id']] = r

    return reactome_node_dict

def getReactomeIdsFromResults(reactome_results):
    """Given the output of getReactomeDataByAcc, extract data meant for building
    graph edges.

    { 'P23219' : [dicts], 'P35354' : [dicts]} => {  'P23219' : [r_id,r_id,...],
                                                    'P35354' : [r_id,r_id,...],
                                                    'ALL'    : [...]}
    """
    acc_reactome_dict = dict()
    all_reactome_ids = set()

    for k1, v1 in reactome_results.items():
        acc_reactome_dict[k1] = list()
        for row in v1:
            all_reactome_ids.add(row['reactome_id'])
            acc_reactome_dict[k1].append(row['reactome_id'])

    acc_reactome_dict['ALL'] = list(all_reactome_ids)
    return acc_reactome_dict

def getChemblVersion():
    data = runQuery("SELECT comments FROM version", None, cursor_factory=psycopg2.extras.DictCursor)
    return data[0]['comments']

def readSources(src_type):
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,src_type+"_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[key.strip()] = val.strip()
    return sources

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

