import psycopg2
import psycopg2.sql as sql
import requests
import json
import os
from psycopg2.extras import NamedTupleCursor
from django.conf import settings

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
    dbhost = 'chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com'
    dbname = 'chembl'
    dbuser = 'chembl'
    dbpass = 'chembl1889'
    
    with psycopg2.connect(host=dbhost, dbname=dbname, user=dbuser, password=dbpass) as conn:
        cur = conn.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
        cur.execute(query, values)
        return cur.fetchall()

def get_chembl_sdfs(chemblIds):
    data = runQuery("""
      SELECT chembl_id||molfile 
      FROM chembl_id_lookup JOIN
           compound_structures ON(entity_id=molregno)
      WHERE chembl_id_lookup.chembl_id IN %s""",(chemblIds,))
    return [row[0] for row in data]

def get_chembl_smiles(chemblIds):
    data = runQuery("""
      SELECT canonical_smiles||' '||chembl_id
      FROM chembl_id_lookup JOIN
           compound_structures ON(entity_id=molregno)
      WHERE chembl_id_lookup.chembl_id IN %s""",(chemblIds,))
    return [row[0] for row in data]

def batchQuery(runQuery,ids, batchSize=1000):
    finalResult = [] 
    print("batching "+str(len(ids))+" ids")
    
    for i in range(0,len(ids),batchSize):
        print("batch "+str(i)+" : "+str(i+batchSize))
        batchResult = runQuery(ids[i:i+batchSize])
        finalResult = finalResult + batchResult
    return finalResult

class MeshIndicationSearch:
    """Class containing MeSH indication search functions and related data"""
    
    mesh_ind_sql = ['chembl_id_lookup.chembl_id',
                    'component_sequences.accession',
                    'drug_mechanism.mechanism_of_action',
                    'drug_mechanism.tid',
                    'component_sequences.component_id',
                    'component_sequences.description',
                    'component_sequences.organism',
                    ]
    
    mesh_ind_list = list()
    for m in mesh_ind_sql:
        m_id = m.replace('.', '__')
        m_sql = m
        mesh_ind_list.append({'id': m_id,
                              'sql': m_sql})
    
    mesh_ind_group_str = 'GROUP BY ' + ', '.join(mesh_ind_sql)
    
    def search(id_type, ids):
        return batchQuery(lambda idBatch: MeshIndicationSearch.__search(id_type,idBatch),ids)

    def __search(id_type, ids):
        """Search for MeSH indications. Basically the original limited-column search.
        
        The idea is to compare the limited columns from the original search with
        corresponding columns in the new AnnotationSearch, and add in the appropriate
        MeSH indication.
        
        Arguments:
        id_type -- String describing type of ID. Must be one of:
                'compound', 'target'
        ids -- Python list of IDs to query
        """
        
        # Adjust query condition based on ID type
        if id_type == 'compound':
            condition = 'chembl_id_lookup.chembl_id = ANY(%s)'
        elif id_type == 'target':
            condition = 'component_sequences.accession = ANY(%s)'
        else:
            raise Exception("MeshIndicationSearch: id_type should be either 'compound' or 'target'")
        
        query = sql.SQL("""SELECT DISTINCT {cols},
                        STRING_AGG(mesh_id||':'||mesh_heading,',') AS mesh_indication
                        FROM chembl_id_lookup
                        JOIN molecule_dictionary USING(chembl_id)
                        JOIN drug_mechanism USING(molregno)
                        JOIN target_components USING(tid)
                        JOIN component_sequences USING(component_id)
                        LEFT JOIN drug_indication USING(molregno)
                        WHERE {condition}
                        {group}
                        ORDER BY 1""")\
                    .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in MeshIndicationSearch.mesh_ind_list),
                    condition=sql.SQL(condition),
                    group=sql.SQL(MeshIndicationSearch.mesh_ind_group_str))
        
        data = runQuery(query, (ids,))
        data = [ d._asdict() for d in data ]
        return data

class DrugIndicationSearch:
    """Class containing drug indication search functions and related data"""
    
    def __init__(self):
        """Initialize a DrugIndicationSearch instance."""
        
        with open(os.path.join(settings.PROJECT_DIR, 'targetsearch/drugind_list.json')) as f:
            self.drugind_list = json.load(f)
    
    def search(self, molregno):
        """Search for drug indications by molregno (compound_structures primary key)"""
        
        query = sql.SQL("""SELECT {cols} FROM drug_indication WHERE molregno = %s""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.drugind_list))
        
        data = runQuery(query, (molregno,))
        data = [ d._asdict() for d in data ]
        return data

class AnnotationSearch:
    """Class containing annotation search functions and related data"""
    
    def __init__(self, id_type):
        """Initialize an AnnotationSearch instance.
        
        AnnotationSearch now needs to be instantiated to allow for per-query
        ordering of columns. The primary purpose of this feature is to set the
        first column to be the groupBy column.
        
        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target'
        """
        
        with open(os.path.join(settings.PROJECT_DIR, 'targetsearch/annotation_list.json')) as f:
            self.annotation_list = json.load(f)
        
        self.id_type = id_type
        
        # Adjust query condition and column ordering based on ID type
        if id_type == 'compound':
            self.condition = 'chembl_id_lookup.chembl_id = ANY(%s)'
            self.groupByIndex = 'annotation__chembl_id_lookup__chembl_id'
        elif id_type == 'target':
            self.condition = 'component_sequences.accession = ANY(%s)'
            self.groupByIndex = 'annotation__component_sequences__accession'
        else:
            raise Exception("AnnotationSearch: id_type should be either 'compound' or 'target'")
        
        # Bring to front the column corresponding to groupByIndex
        for col in self.annotation_list:
            if col['id'] == self.groupByIndex:
                self.annotation_list.remove(col)
                self.annotation_list.insert(0, col)
                break
    
    def search(self, ids):
        return batchQuery(lambda idBatch: self.__search(idBatch),ids)
    def __search(self, ids):
        """Search for targets by annotations, returning all database columns.
        
        By returning all columns, we allow the user to decide which columns are
        desired, by means of the Django/DataTables interface.
        
        Do note that the argument style is different than byAnnotations().
        
        Arguments:
        ids -- Python list of IDs to query
        """
        
        query = sql.SQL("""SELECT DISTINCT {cols} FROM chembl_id_lookup
                        JOIN molecule_dictionary USING(chembl_id)
                        JOIN drug_mechanism USING(molregno)
                        JOIN target_components USING(tid)
                        JOIN component_sequences USING(component_id)
                        WHERE {condition}
                        ORDER BY 1""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.annotation_list),
                        condition=sql.SQL(self.condition))
        
        data = runQuery(query, (ids,))
        data = [ d._asdict() for d in data ]
        return data
    
    def search_grouped(self, ids):
        """Search for targets by annotations, returning all database columns.
        
        search_grouped() performs everything in search(), plus the additional
        step of grouping the rows with groupBy(). This is closer to the original
        byAnnotations() output structure. However, the original search()
        structure may be easier to work with in certain scenarios.
        
        Arguments:
        ids -- Python list of IDs to query"""
        
        data = self.search(ids)
        return groupBy(lambda t: t[self.groupByIndex], data)

class AnnotationWithMeshSearch:
    """Special class for adding MeSH indication data to annotation search"""
    
    def __init__(self, id_type):
        """Initialize an AnnotationWithMeshSearch instance.
        
        AnnotationWithMeshSearch needs to be instantiated to allow for
        per-query ordering of columns. The primary purpose of this feature is
        to set the first column to be the groupBy column.
        
        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target'
        """
        
        self.id_type = id_type
        self.myAnnotationSearch = AnnotationSearch(id_type)
        self.annotation_list = self.myAnnotationSearch.annotation_list.copy()
        self.annotation_list.append({'id': 'annotation__mesh_indication',
                                     'sql': "string_agg(mesh_id||':'||mesh_heading,',')",
                                     'name': 'MeSH Indication',
                                     'desc': "Aggregated list of MeSH identifiers and heading terms",
                                     'visible': False,
                                     'url': None,
                                     'childrow': True,
                                     })
    
    def search(self, ids):
        return batchQuery(lambda idBatch: self.__search(idBatch),ids)
    def __search(self, ids):
        annotation_data = self.myAnnotationSearch.search(ids)
        mesh_ind_data = MeshIndicationSearch.search(self.id_type, ids)
        
        for a in annotation_data:
            mesh_ind = None
            for m in mesh_ind_data:
                if  a['annotation__chembl_id_lookup__chembl_id'] == m['chembl_id_lookup__chembl_id'] and \
                    a['annotation__component_sequences__accession'] == m['component_sequences__accession'] and \
                    a['annotation__drug_mechanism__mechanism_of_action'] == m['drug_mechanism__mechanism_of_action'] and \
                    a['annotation__drug_mechanism__tid'] == m['drug_mechanism__tid'] and \
                    a['annotation__component_sequences__component_id'] == m['component_sequences__component_id'] and \
                    a['annotation__component_sequences__description'] == m['component_sequences__description'] and \
                    a['annotation__component_sequences__organism'] == m['component_sequences__organism']:
                        mesh_ind = m['mesh_indication']
                        break
            a['annotation__mesh_indication'] = mesh_ind
        
        return annotation_data
    
    def search_grouped(self, ids):
        data = self.search(ids)
        return groupBy(lambda t: t[self.myAnnotationSearch.groupByIndex], data)

class ActivitySearch:
    """Class containing activity search functions and related data"""
    
    def __init__(self, id_type):
        """Initialize an ActivitySearch instance.
        
        ActivitySearch now needs to be instantiated to allow for per-query
        ordering of columns. The primary purpose of this feature is to set the
        first column to be the groupBy column.
        
        Arguments:
        id_type -- String describing type of ID. Must be one of:
                   'compound', 'target'
        """
        
        with open(os.path.join(settings.PROJECT_DIR, 'targetsearch/activity_list.json')) as f:
            self.activity_list = json.load(f)
        
        self.id_type = id_type
        
        # Adjust query condition and column ordering based on ID type
        if id_type == 'compound':
            self.condition = 'chembl_id_lookup.chembl_id = ANY(%s)'
            self.groupByIndex = 'activity__chembl_id_lookup__chembl_id'
        elif id_type == 'target':
            self.condition = "chembl_id_lookup.entity_type='COMPOUND' AND component_sequences.accession = ANY(%s)"
            self.groupByIndex = 'activity__component_sequences__accession'
        else:
            raise Exception("ActivitySearch: id_type should be either 'compound' or 'target'")
        
        # Bring to front the column corresponding to groupByIndex
        for col in self.activity_list:
            if col['id'] == self.groupByIndex:
                self.activity_list.remove(col)
                self.activity_list.insert(0, col)
                break
    
    def search(self, ids):
        return batchQuery(lambda idBatch: self.__search(idBatch),ids)
    def __search(self, ids):
        """Search for targets by activity, returning all database columns.
        
        By returning all columns, we allow the user to decide which columns are
        desired, by means of the Django/DataTables interface.
        
        Do note that the argument style is different than byActivity().
        
        Arguments:
        ids -- Python list of IDs to query
        """
        
        query = sql.SQL("""SELECT {cols} FROM chembl_id_lookup
                        JOIN activities ON(entity_id = molregno)
                        JOIN molecule_dictionary USING(molregno)
                        JOIN assays USING(assay_id)
                        JOIN target_components USING(tid)
                        JOIN component_sequences USING(component_id)
                        WHERE {condition}
                        ORDER BY 1""")\
                .format(cols=sql.SQL(', ').join(sql.SQL("{} AS {}".format(c['sql'], c['id'])) for c in self.activity_list),
                        condition=sql.SQL(self.condition))
        
        data = runQuery(query, (ids,))
        data = [ d._asdict() for d in data ]
        return data
    
    def search_grouped(self, ids):
        data = self.search(ids)
        return groupBy(lambda t: t[self.groupByIndex], data)

def mapToChembl(unknownIds, sourceId):
    chemblIds = set()
    for unknownId in unknownIds:
        req = requests.get("https://www.ebi.ac.uk/unichem/rest/src_compound_id/"+unknownId+"/"+str(sourceId)+"/1")
        #print("\n\nreq.json: "+str(req.json()))
        result = req.json()
        #print("got result: "+str(result))
        if "error" in result: 
            raise Exception(result["error"])
        if isinstance(result,list) and len(result) > 0 and ("src_compound_id" in result[0]) :
            chemblIds.add(result[0]["src_compound_id"])

    #print("final chembl ids: "+str(chemblIds))

    return tuple(chemblIds)

def mapToUniprot(unknownIds, sourceId):
    url = "https://www.uniprot.org/uploadlists/"
    params = {
            "from": sourceId,
            "to": "ACC",
            "format": "tab",
            "query": " ".join(unknownIds)
            }

    results = requests.get(url,params=params).text
    return tuple([ line.split("\t")[1]  for line in results.splitlines() if  "From\tTo" not in line])



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
    sqlQuery = sql.SQL("""
               SELECT DISTINCT chembl_id, synonyms 
               FROM molecule_synonyms 
                JOIN molecule_dictionary USING(molregno) 
                WHERE synonyms ilike %s """)
	
    return runQuery(sqlQuery,("%"+nameQuery+"%",))


#getUniChemSources()
#mapToChembl(['DB00829','DB00945'],2)

#print(chemblTargetAccessionsByActivity(('CHEMBL26',)))
#print(chemblTargetAccessionsByAnnotations(('CHEMBL25',)))
