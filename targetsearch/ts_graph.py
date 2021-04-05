import json
from .helpers import ActivitySearch, AnnotationWithDrugIndSearch

def pj(obj):
    print(json.dumps(obj, indent=2))

def get_anno_data(id_type, ids):
    myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, ids)
    anno_data = list()
    #for r in myAnnotationSearch.query_data_dict:
    for r in myAnnotationSearch.get_results():
        d = dict()
        #d.update(r)
        d['chembl_id'] = r['annotation__chembl_id_lookup__chembl_id']
        d['entity_id'] = r['annotation__chembl_id_lookup__entity_id']
        d['status'] = r['annotation__chembl_id_lookup__status']
        d['pref_name'] = r['annotation__molecule_dictionary__pref_name']
        d['mech_action'] = r['annotation__drug_mechanism__mechanism_of_action']
        d['tid'] = r['annotation__drug_mechanism__tid']
        d['component_id'] = r['annotation__component_sequences__component_id']
        d['acc_id'] = r['annotation__component_sequences__accession']
        d['desc'] = r['annotation__component_sequences__description']
        d['org'] = r['annotation__component_sequences__organism']
        anno_data.append(d)
    return anno_data

def get_act_data(id_type, ids):
    myActivitySearch = ActivitySearch(id_type, ids)
    act_data = list()
    for r in myActivitySearch.get_results():
        d = dict()
        d['chembl_id'] = r['activity__chembl_id_lookup__chembl_id']
        d['entity_id'] = r['activity__chembl_id_lookup__entity_id']
        d['std_rel'] = r['activity__activities__standard_relation']
        d['std_val'] = r['activity__activities__standard_value']
        d['std_units'] = r['activity__activities__standard_units']
        d['std_flag'] = r['activity__activities__standard_flag']
        d['std_type'] = r['activity__activities__standard_type']
        d['act_comment'] = r['activity__activities__activity_comment']
        d['dv_comment'] = r['activity__activities__data_validity_comment']
        d['pot_dup'] = r['activity__activities__potential_duplicate']
        d['pchembl_val'] = r['activity__activities__pchembl_value']
        d['upper_val'] = r['activity__activities__upper_value']
        d['std_upper_val'] = r['activity__activities__standard_upper_value']
        d['pref_name'] = r['activity__molecule_dictionary__pref_name']
        d['assay_id'] = r['activity__assays__assay_id']
        d['assay_desc'] = r['activity__assays__description']
        d['assay_chembl_id'] = r['activity__assays__chembl_id']
        d['component_type'] = r['activity__component_sequences__component_type']
        d['acc_id'] = r['activity__component_sequences__accession']
        d['desc'] = r['activity__component_sequences__description']
        d['org'] = r['activity__component_sequences__organism']
        act_data.append(d)
    return act_data

def make_anno_graph(query_data, table_name = 'annotation'):
    node_id_next = 1
    nodes = list()
    edges = list()
    compounds = dict()
    targets = dict()
    cmp_nid_lookup = dict() # d[node_id] -> chembl_id, d[chembl_id] -> node_id
    tgt_nid_lookup = dict() # d[node_id] -> acc, d[acc] -> node_id

    # Gather unique compound info
    for r in query_data:
        chembl_id = r['chembl_id']

        if chembl_id not in compounds:
            cmp_info = dict()
            for k in ['chembl_id', 'pref_name']:
                cmp_info[k] = r[k]
            compounds[chembl_id] = cmp_info

    # Gather unique target info
    for r in query_data:
        acc_id = r['acc_id']

        if acc_id not in targets:
            tgt_info = dict()
            if table_name == 'annotation':
                keys = ['acc_id', 'desc', 'org']
            elif table_name == 'activity':
                keys = ['assay_id', 'std_type', 'std_rel', 'std_val', 'pref_name',
                        'component_type', 'acc_id', 'desc', 'org']
            else:
                raise Exception('table_name should be annotation or activity')
            for k in keys:
                tgt_info[k] = r[k]
            targets[acc_id] = tgt_info

    # Add compounds to node list
    for chembl_id, cmp_info in compounds.items():
        node_id = node_id_next
        node_id_next += 1

        n = dict()
        n['id'] = node_id
        n['label'] = chembl_id
        n['shape'] = 'box'
        n['type'] = 'compound'
        n['group'] = 1
        #n['color'] = '#abcdef'

        t = '{chembl_id}\n{pref_name}'
        n['title'] = t.format(**cmp_info)

        nodes.append(n)
        cmp_nid_lookup[node_id] = chembl_id
        cmp_nid_lookup[chembl_id] = node_id

    # Add annotation targets to node list
    for acc_id, tgt_info in targets.items():
        node_id = node_id_next
        node_id_next += 1

        n = dict()
        n['id'] = node_id
        n['label'] = acc_id
        n['shape'] = 'ellipse'
        n['type'] = 'target'
        n['group'] = 2
        #n['color'] = '#abcdef'

        t = '{acc_id}\n{desc}\n{org}'
        n['title'] = t.format(**tgt_info)

        nodes.append(n)
        tgt_nid_lookup[node_id] = acc_id
        tgt_nid_lookup[acc_id] = node_id

    # Create edges
    for r in query_data:
        cmp_node = cmp_nid_lookup[r['chembl_id']]
        tgt_node = tgt_nid_lookup[r['acc_id']]

        e = dict()
        e['from'] = cmp_node
        e['to'] = tgt_node
        edges.append(e)

    return (nodes, edges)

def make_act_graph(query_data):
    node_id_next = 1
    nodes = list()
    edges = list()
    compounds = dict()
    targets = dict()
    cmp_nid_lookup = dict() # d[node_id] -> chembl_id, d[chembl_id] -> node_id
    tgt_nid_lookup = dict() # d[node_id] -> acc, d[acc] -> node_id
    assays_dict = dict()

    # Find unique compounds
    for r in query_data:
        chembl_id = r['chembl_id']

        if chembl_id not in compounds:
            cmp_info = dict()
            for k in ['chembl_id', 'pref_name']:
                cmp_info[k] = r[k]
            compounds[chembl_id] = cmp_info

    # Gather activity target and assays info
    for r in query_data:
        acc_id = r['acc_id']

        if acc_id not in targets:
            tgt_info = dict()
            tgt_keys = ['component_type', 'acc_id', 'desc', 'org']
            for k in tgt_keys:
                tgt_info[k] = r[k]
            targets[acc_id] = tgt_info

        assay_info = dict()
        assay_keys = ['assay_chembl_id', 'std_type', 'std_rel', 'std_val'] #'component_type',
        #        'acc_id', 'desc', 'org']
        for k in assay_keys:
            #tgt_info_list[k] = r[k]
            assay_info[k] = r[k]

        # An accession ID may have multiple assays associated with it.
        #if acc_id not in targets:
        #    targets[acc_id] = list()
        #targets[acc_id].append(tgt_info_list)
        if acc_id not in assays_dict:
            assays_dict[acc_id] = list()
        assays_dict[acc_id].append(assay_info)


    # Add compounds to node list
    for chembl_id, cmp_info in compounds.items():
        node_id = node_id_next
        node_id_next += 1

        n = dict()
        n['id'] = node_id
        n['label'] = chembl_id
        n['shape'] = 'box'
        n['type'] = 'compound'
        n['group'] = 1
        #n['color'] = '#abcdef'

        t = '{chembl_id}\n{pref_name}'
        n['title'] = t.format(**cmp_info)

        nodes.append(n)
        cmp_nid_lookup[node_id] = chembl_id
        cmp_nid_lookup[chembl_id] = node_id

    # Add activity targets to node list
    for acc_id, tgt_info in targets.items():
        node_id = node_id_next
        node_id_next += 1

        n = dict()
        n['id'] = node_id
        n['label'] = acc_id
        n['shape'] = 'ellipse'
        n['type'] = 'target'
        n['group'] = 2
        #n['color'] = '#abcdef'

        title_list = list()
        t = '{component_type}\n{acc_id}\n{desc}\n{org}'
        title_list.append(t.format(**tgt_info))
        #for tgt_info in tgt_info_list:
            #t = 'Assay ID: {assay_id}\n' +\
            #'{std_type} {std_rel} {std_val}\n\n' +\
            #'{component_type}\n{acc_id}\n{desc}\n{org}'
        for assay_info in assays_dict[acc_id]:
            t = 'Assay: {assay_chembl_id}\n{std_type} {std_rel} {std_val}'
            title_list.append(t.format(**assay_info))
        n['title'] = '\n----------\n'.join(title_list)

        nodes.append(n)
        tgt_nid_lookup[node_id] = acc_id
        tgt_nid_lookup[acc_id] = node_id

    # Create edges
    for r in query_data:
        cmp_node = cmp_nid_lookup[r['chembl_id']]
        tgt_node = tgt_nid_lookup[r['acc_id']]

        e = dict()
        e['from'] = cmp_node
        e['to'] = tgt_node
        edges.append(e)

    return (nodes, edges)
