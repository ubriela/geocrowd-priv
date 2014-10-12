# To change this template, choose Tools | Templates
# and open the template in the editor.

import unittest

from geocrowd.HT_pure import HT_pure


class HTree_TestCase(unittest.TestCase):
    # def setUp(self):
    #    self.foo = HTree_()
    #

    #def tearDown(self):
    #    self.foo.dispose()
    #    self.foo = None

    def test_hTree_(self):
        #assert x != y;
        ht_pure = HT_pure()
        height = 2
        array = [1, 2, 3, 4, 5, 6, 7]
        left, right = 0, 9
        budget_s = [0 for _ in range(3)]
        x = ht_pure.recursiveMedians(height, array, left, right, budget_s)
        print x


if __name__ == '__main__':
    unittest.test_hTree_

