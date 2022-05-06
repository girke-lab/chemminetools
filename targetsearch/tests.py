from django.test import TestCase, SimpleTestCase
from django.urls import reverse

import targetsearch.helpers as helpers
import targetsearch.views as views

# Create your tests here.

class TsFilterTestCase(SimpleTestCase):
    def test_getGoIdsByAccFlat(self):
        root_nodes = {'GO:0008150', 'GO:0005575', 'GO:0003674'} # process, component, function

        gf = helpers.getGoIdsByAccFlat(['P23219'], remove_roots=True, add_missing_links=True)
        a = set(gf['ALL'])
        b = set(gf['P23219'])
        self.assertEqual(a, b, "missing intermediate nodes in flat mode")

        # Expected values for P02918
        # NOTE: This may need to be updated with newer versions of ChEMBL
        correct_acc_go_edges = {'GO:0016787', 'GO:0071554', 'GO:0006508', 'GO:0016740', 'GO:0009636', 'GO:0046677', 'GO:0016020', 'GO:0009058', 'GO:0005515', 'GO:0043167', 'GO:0003824', 'GO:0005886', 'GO:0008233', 'GO:0008152', 'GO:0016757', 'GO:0016021'}
        test_acc_go_edges = helpers.getGoIdsByAccFlat(['P02918'], remove_roots=True, add_missing_links=True)
        for rn in root_nodes:
            #self.assertNotIn(rn, test_acc_go_edges['P02918'], f"root node in acc_go_edges: {rn}")
            self.assertNotIn(rn, test_acc_go_edges['P02918'], "root node in acc_go_edges: {}".format(rn))
        self.assertEqual(len(correct_acc_go_edges), len(test_acc_go_edges['P02918']), "acc_go_edges has incorrect number of items")
        test_acc_go_edges_set = set(test_acc_go_edges['P02918'])
        self.assertEqual(len(test_acc_go_edges['P02918']), len(test_acc_go_edges_set), "acc_go_edges contains duplicates")
        self.assertEqual(correct_acc_go_edges, test_acc_go_edges_set, "acc_go_edges contents are incorrect")

    def test_tsAnnoFilter1(self):
        t = views.tsAnnoFilter1_debug(['CHEMBL174'])

        correct_age_keys = {'P0AD65', 'P0AEB2', 'P02918', 'P08506', 'P0AD68', 'P24228', 'P02919', 'ALL'}
        self.assertEqual(correct_age_keys, set(t['acc_go_edges'].keys()), "incorrect set of keys in acc_go_edges")

class TargetSearchViewTests(TestCase):
    def test_ts_main(self):
        response = self.client.get(reverse('newTS'))
        self.assertEqual(response.status_code, 200)

    def test_tsAnnoFilter1_blank(self):
        response = self.client.get(reverse('tsAnnoFilter1'))
        self.assertEqual(response.status_code, 200)

    def test_tsAnnoFilter1_with_ids(self):
        response = self.client.get(reverse('tsAnnoFilter1'), {'ids': ['CHEMBL25', 'CHEMBL112', 'CHEMBL154', 'CHEMBL521']})
        self.assertEqual(response.status_code, 200)

    def test_tsAnnoFilter1_no_results(self):
        response = self.client.get(reverse('tsAnnoFilter1'), {'ids': ['CHEMBL0']})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Search returned no annotation data")
