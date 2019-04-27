import psycopg2

def tupleArray2Dict(tuples):
    dict = {}
    for tuple in tuples:
        a = dict.get(tuple[0],[])
        a.append(tuple[1])
        dict[tuple[0]] = a
    return dict

def runQuery(query,values):
    conn = psycopg2.connect(host = "chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com",
            user="chembl",dbname="chembl",password="chembl1889")
    cur = conn.cursor()
    cur.execute(query,values)
    return cur.fetchall()

def chemblTargetAccessions(chemblIds):
    data = runQuery("""
            select distinct chembl_id_lookup.chembl_id, accession
                from chembl_id_lookup
                        join activities on(entity_id = molregno)
                        join assays using(assay_id)
                        join target_components using(tid)
                        join component_sequences using(component_id)
                where 
                        chembl_id_lookup.chembl_id IN ( %s )
                order by 1
            """,(chemblIds,))
    return tupleArray2Dict(data)

def accessionToChembl(accessionIds):
    data = runQuery("""
            select distinct accession, chembl_id_lookup.chembl_id
                from chembl_id_lookup
                        join activities on(entity_id = molregno)
                        join assays using(assay_id)
                        join target_components using(tid)
                        join component_sequences using(component_id)
                where
                        accession IN ( %s )
                order by 1
            """,[accessionIds])
    return tupleArray2Dict(data)

def mapToChembl(unknownIds):
    chemblIds = set()
    for unknownId in unknownIds:
        req = requests.get("https://www.ebi.ac.uk/unichem/rest/orphanIdMap/"+unknownId+"/1")
        data = [ source[0]["src_compound_id"]  for source in req.json().values()]
        #print("got result: "+str(data))
        chemblIds |= set(data)

    #print("final chembl ids: "+str(chemblIds))

    return list(chemblIds)
