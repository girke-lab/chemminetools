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
from psycopg2.extras import NamedTupleCursor
from django.conf import settings

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

def runQuery(query, values):
    dbhost = settings.CHEMBL_DB['DBHOST']
    dbname = settings.CHEMBL_DB['DBNAME']
    dbuser = settings.CHEMBL_DB['DBUSER']
    dbpass = settings.CHEMBL_DB['DBPASS']

    with psycopg2.connect(host=dbhost, dbname=dbname, user=dbuser, password=dbpass) as conn:
        cur = conn.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
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

class SearchBase:
    """Base class for specialized search sub-classes. Should not be used directly."""

    def __init__(self):
        # Dict-based column data. Use load_json_info or manually define/edit
        self.table_info = None

        # Raw query results returned from fetchall()
        self.query_data = None

        # Name of column by which to group results
        self.groupByIndex = None

    def load_json_info(self, filename):
        """Method for loading table_info from a JSON file. Should be called
        inside __init__ in subclasses."""
        self.table_info = SearchBase._load_json_info(filename)

    def _load_json_info(filename):
        """Static function variant of load_json_info. Use if for some reason
        you don't want to immediately assign table_info."""
        with open(os.path.join(settings.PROJECT_DIR, filename)) as f:
            return json.load(f)

    def run_query(self, query, values):
        """Method for running a general purpose query. Stores results in
        query_data. Should be called inside __init__ in subclasses.

        Does not handle large argument sets. For that, either use run_batch_query
        or wrap this method inside batchQuery."""
        self.query_data = SearchBase._run_query(query, values)

    def _run_query(query, values):
        """Static function variant of run_query. Use if for some reason you
        don't want to immediately assign query_data."""
        return runQuery(query, values)

    def run_batch_query(self, query, ids):
        """Convenience method for running a batch query. Stores results in
        query_data. Should be called inside __init__ in subclasses.

        Only applicable when there is one SQL field to fill."""
        self.query_data = SearchBase._run_batch_query(query, ids)

    def _run_batch_query(query, ids):
        """Static function variant of run_batch_query. Use if for some reason
        you don't want to immediately assign query_data."""
        return batchQuery(lambda idBatch: runQuery(query, (idBatch,)), ids)

    def get_results(self):
        """Return query data as an OrderedDict list."""
        return [ d._asdict() for d in self.query_data ]

    def get_grouped_results(self):
        """Return query data in a grouped structure."""
        return groupBy(lambda t: t[self.groupByIndex], self.get_results())

class DrugIndicationSearch(SearchBase):
    """Class containing drug indication search functions and related data"""

    def __init__(self, molregno):
        """Initialize a DrugIndicationSearch instance.

        Arguments:
        molregno -- String for searching drug_indication by molregno.
        """

        super()
        self.molregno = molregno

        # Load column information from JSON
        super().load_json_info('targetsearch/drugind_info.json')

        # Run the query
        query = sql.SQL("""SELECT {cols} FROM drug_indication WHERE molregno = %s""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.table_info))
        super().run_query(query, (molregno,))

class AnnotationSearch(SearchBase):
    """Class containing annotation search functions and related data"""

    def __init__(self, id_type, ids):
        """Initialize an AnnotationSearch instance.

        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target', 'homolog-target'
        ids -- Python list of IDs to query
        """

        super()

        # Load column information from JSON
        super().load_json_info('targetsearch/annotation_info.json')

        self.id_type = id_type

        # Adjust query condition and column ordering based on ID type
        if id_type == 'compound':
            self.condition = 'chembl_id_lookup.chembl_id = ANY(%s)'
            self.groupByIndex = 'annotation__chembl_id_lookup__chembl_id'
        elif id_type in ['target', 'homolog-target']:
            self.condition = 'component_sequences.accession = ANY(%s)'
            self.groupByIndex = 'annotation__component_sequences__accession'
        else:
            raise Exception("AnnotationSearch: id_type should be either 'compound', 'target', or 'homolog-target'")

        # Bring to front the column corresponding to groupByIndex
        for col in self.table_info:
            if col['id'] == self.groupByIndex:
                self.table_info.remove(col)
                self.table_info.insert(0, col)
                break

        # Run the query
        query = sql.SQL("""SELECT DISTINCT {cols} FROM chembl_id_lookup
                        JOIN molecule_dictionary USING(chembl_id)
                        JOIN drug_mechanism USING(molregno)
                        JOIN target_components USING(tid)
                        JOIN component_sequences USING(component_id)
                        WHERE {condition}
                        ORDER BY 1""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.table_info),
                        condition=sql.SQL(self.condition))
        super().run_batch_query(query, ids)

        # Extract ChEMBL ID and molregno from query_data
        temp_molregno_set = set()
        for r in self.query_data:
            pair = (r.annotation__chembl_id_lookup__chembl_id, r.annotation__molecule_dictionary__molregno)
            temp_molregno_set.add(pair)

        self.molregno_set = set()
        self.molregno_to_chembl = dict()
        for chembl, molregno in temp_molregno_set:
            self.molregno_set.add(molregno)
            self.molregno_to_chembl[molregno] = chembl

class AnnotationWithDrugIndSearch(AnnotationSearch):
    """Special class for adding drug indication data to Annotation Search"""

    def __init__(self, id_type, ids):
        """Initialize an AnnotationWithDrugIndSearch instance.

        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target', 'homolog-target'
        ids -- Python list of IDs to query
        """

        super().__init__(id_type, ids)

        # Save a copy of the AnnotationSearch results as an OrderedDict.
        # This is the authoritative data source for objects of this class.
        self.query_data_dict = super().get_results()

        # Add extended table information
        new_table_info = [
            {
            'id': 'annotation__drugind_summary__mesh_indications',
            'sql': None,
            'table': 'Drug Indication Summary',
            'name': 'MeSH Indications',
            'desc': 'Concatenated list of MeSH indications.',
            'visible': False,
            },
            {
            'id': 'annotation__drugind_summary__efo_indications',
            'sql': None,
            'table': 'Drug Indication Summary',
            'name': 'EFO Indications',
            'desc': 'Concatenated list of EFO indications.',
            'visible': False,
            },
            {
            'id': 'annotation__drugind_summary__show_drugind_table',
            'sql': None,
            'table': 'Drug Indication Summary',
            'name': 'Drug Indication Table',
            'desc': 'Button to view Drug Indication table modal.',
            'visible': True,
            'export': False,
            'html': '<button type="button" class="btn btn-primary" onclick="showDrugIndModal(\'{}\')">View</button>',
            }]
        self.table_info += new_table_info

        # Add blank cells for the new table columns
        for row in self.query_data_dict:
            for n in new_table_info:
                row[n['id']] = None

        # Prepare Drug Indication search data
        self.drugind_objs = dict()
        for m in self.molregno_set:
            c = self.molregno_to_chembl[m]
            self.drugind_objs[c] = DrugIndicationSearch(m)
            self.add_drugind_data(self.drugind_objs[c])

    def add_drugind_data(self, drugind_obj):
        """Helper method. Extract data from a prepared DrugIndicationSearch
        object and add to this object's Annotation Search table. Can be run
        multiple times, once for each molregno."""

        # Get results from DrugIndicationSearch
        drugind_data_dict = drugind_obj.get_results()

        mesh_indications = list()
        efo_indications = list()

        # Prepare MeSH and EFO indication substrings
        for d in drugind_data_dict:
            mesh_id = d['drug_indication__mesh_id']
            mesh_heading = d['drug_indication__mesh_heading']
            efo_id = d['drug_indication__efo_id']
            efo_term = d['drug_indication__efo_term']

            if mesh_id != None and mesh_heading != None:
                mesh_indications.append( "[{}] {}".format(mesh_id, mesh_heading) )
            if efo_id != None and efo_term != None:
                efo_indications.append( "[{}] {}".format(efo_id, efo_term) )

        # Insert data into main table
        for row in self.query_data_dict:
            if row['annotation__molecule_dictionary__molregno'] == drugind_obj.molregno:
                row['annotation__drugind_summary__mesh_indications'] = "; ".join(mesh_indications)
                row['annotation__drugind_summary__efo_indications'] = "; ".join(efo_indications)
                row['annotation__drugind_summary__show_drugind_table'] = self.molregno_to_chembl[drugind_obj.molregno]

    def get_results(self):
        """Return query data as an OrderedDict list."""
        return self.query_data_dict

class ActivitySearch(SearchBase):
    """Class containing activity search functions and related data"""

    def __init__(self, id_type, ids):
        """Initialize an ActivitySearch instance.

        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target', 'homolog-target'
        ids -- Python list of IDs to query
        """

        super()

        # Load column information from JSON
        super().load_json_info('targetsearch/activity_info.json')

        self.id_type = id_type

        # Adjust query condition and column ordering based on ID type
        if id_type == 'compound':
            self.condition = 'chembl_id_lookup.chembl_id = ANY(%s)'
            self.groupByIndex = 'activity__chembl_id_lookup__chembl_id'
        elif id_type in ['target', 'homolog-target']:
            self.condition = "chembl_id_lookup.entity_type='COMPOUND' AND component_sequences.accession = ANY(%s)"
            self.groupByIndex = 'activity__component_sequences__accession'
        else:
            raise Exception("ActivitySearch: id_type should be either 'compound', 'target', or 'homolog-target'")

        # Bring to front the column corresponding to groupByIndex
        for col in self.table_info:
            if col['id'] == self.groupByIndex:
                self.table_info.remove(col)
                self.table_info.insert(0, col)
                break

        # Run the query
        query = sql.SQL("""SELECT {cols} FROM chembl_id_lookup
                        JOIN activities ON(entity_id = molregno)
                        JOIN molecule_dictionary USING(molregno)
                        JOIN assays USING(assay_id)
                        JOIN target_components USING(tid)
                        JOIN component_sequences USING(component_id)
                        WHERE {condition}
                        ORDER BY 1""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.table_info),
                        condition=sql.SQL(self.condition))
        super().run_batch_query(query, ids)

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





                

