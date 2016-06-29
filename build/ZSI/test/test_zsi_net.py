#!/usr/bin/env python
import unittest
import test_t1
import test_t2
import test_t3
import test_t4
import test_t5
import test_t6
import test_t7
import test_t8
import test_t9

def makeTestSuite():
    suite1 = test_t1.makeTestSuite()
    suite2 = test_t2.makeTestSuite()
    suite3 = test_t3.makeTestSuite()
    suite4 = test_t4.makeTestSuite()
    suite5 = test_t5.makeTestSuite()
    suite6 = test_t6.makeTestSuite()
    suite7 = test_t7.makeTestSuite()
    suite8 = test_t8.makeTestSuite()
    suite9 = test_t9.makeTestSuite()
    t = (suite1, suite2, suite3, suite4, suite5, suite6, suite7, suite8, suite9)
    suite = unittest.TestSuite(t)
    return suite
def main():
    unittest.main(defaultTest="makeTestSuite")
    suite = unittest.TestSuite()

if __name__ == "__main__" : main()
