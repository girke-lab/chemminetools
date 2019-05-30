import psycopg2
import requests
from psycopg2.extras import NamedTupleCursor


def tupleArray2Dict(tuples):
    dict = {}
    for tuple in tuples:
        a = dict.get(tuple[0],[])
        a.append(tuple[1])
        dict[tuple[0]] = a
    return dict

def groupBy(keyFn, tuples):
    dict = {}
    for tuple in tuples:
        key = keyFn(tuple)
        a = dict.get(key,[])
        a.append(tuple)
        dict[key] = a

    return dict

def runQuery(query,values):
    conn = psycopg2.connect(host = "chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com",
            user="chembl",dbname="chembl",password="chembl1889")
    cur = conn.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
    cur.execute(query,values)
    return cur.fetchall()

def chemblTargetAccessions(chemblIds):
    data = runQuery("""
            select distinct chembl_id_lookup.chembl_id, accession,
                    activities.molregno, pref_name, activity_id, assays.chembl_id AS chembl_assay_id, 
                    component_sequences.description, organism, 
                    activities.standard_value, activities.standard_units, 
                    activities.standard_flag, activities.standard_type
                from chembl_id_lookup
                        join activities on(entity_id = molregno)
                        join molecule_dictionary using(molregno)
                        join assays using(assay_id)
                        join target_components using(tid)
                        join component_sequences using(component_id)
                where 
                        chembl_id_lookup.chembl_id in %s
                order by 1
            """,(chemblIds,))
    #return tupleArray2Dict(data)
    return groupBy(lambda t: t[0], data)
def chemblTargetAccessionsByAnnotations(chemblIds):
    data = runQuery("""
            select distinct chembl_id,  accession,
                    mechanism_of_action, tid,component_id, description, 
                    organism, mesh_id, mesh_heading
                from chembl_id_lookup
                        join molecule_dictionary using(chembl_id)
                        join drug_mechanism using(molregno)
                        join drug_indication using(molregno)
                        join target_components using(tid)
                        join component_sequences using(component_id)
                where 
                        chembl_id_lookup.chembl_id in %s
                order by 1
            """,(chemblIds,))
    #return tupleArray2Dict(data)
    return groupBy(lambda t: t[0], data)


def accessionToChembl(accessionIds):
    data = runQuery("""
            select distinct accession, chembl_id_lookup.chembl_id
                from chembl_id_lookup
                        join activities on(entity_id = molregno)
                        join assays using(assay_id)
                        join target_components using(tid)
                        join component_sequences using(component_id)
                where
                        accession in %s
                order by 1
            """,(accessionIds,))
    #return tupleArray2Dict(data)
    return groupBy(lambda t: t[0], data)

def mapToChembl(unknownIds, sourceId):
    chemblIds = set()
    for unknownId in unknownIds:
        req = requests.get("https://www.ebi.ac.uk/unichem/rest/src_compound_id/"+unknownId+"/"+str(sourceId)+"/1")
        #print("\n\nreq.json: "+str(req.json()))
        result = req.json()
        #print("got result: "+str(result))
        if len(result) > 0 and ("src_compound_id" in result[0]) :
            chemblIds.add(result[0]["src_compound_id"])

    #print("final chembl ids: "+str(chemblIds))

    return tuple(chemblIds)

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


#getUniChemSources()
#mapToChembl(['DB00829','DB00945'],2)

#print(chemblTargetAccessions(('CHEMBL26',)))
#print(chemblTargetAccessionsByAnnotations(('CHEMBL25',)))