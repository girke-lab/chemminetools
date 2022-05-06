# chembl_search.py - ChEMBL search class definitions

from .helpers import *

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

        super().__init__()
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

        super().__init__()

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
            },
            {
            'id': 'annotation__extanno__show_extanno',
            'sql': None,
            'table': 'Extended Annotations',
            'name': 'Extended Annotations',
            'desc': 'Link to view extended annotations.',
            'visible': False,
            'export': False,
            #'html': '<a class="btn btn-primary" role="button" target="_blank" href="/targetsearch/extannobychembl/{}/">View</a>',
            'extanno': True,
            }
            ]
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

        # Prepare Extended Annotation data

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
                #row['annotation__extanno__show_extanno'] = self.molregno_to_chembl[drugind_obj.molregno]
                extanno_dict = dict()
                chembl_id = self.molregno_to_chembl[drugind_obj.molregno]
                extanno_dict['chembl_id'] = chembl_id
                extanno_check = checkExtAnno(chembl_id)
                if extanno_check is None:
                    extanno_dict['extannos'] = None
                else:
                    extanno_dict['extannos'] = set()
                    for k, v in extanno_check.items():
                        if v is not None:
                            extanno_dict['extannos'].add(k)
                row['annotation__extanno__show_extanno'] = extanno_dict

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

        super().__init__()

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
