from django.test import TestCase
import unittest # for test cases that do not need the Django test database

import targetsearch.helpers as helpers

# Create your tests here.

class TsFilterTestCase(unittest.TestCase):
    def test_getGoIdsByAcc(self):
        gf = helpers.getGoIdsByAcc(['P23219'], flat=True)
        a = set(gf['ALL'])
        b = set(gf['P23219'])
        self.assertEqual(a, b, "missing intermediate nodes in flat mode")
