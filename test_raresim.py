import unittest
from rareSim import sparse
from header import *
import random

class testRaresim(unittest.TestCase):
    
    def test_readExpected(self):
        bins = read_expected('./testData/testBins.txt')
        self.assertEqual(bins[2], (3,5,2))

    def test_getBin(self):
        bins = read_expected('./testData/testBins.txt')
        self.assertEqual(get_bin(bins, 4), 2)

    def test_rowNum(self):
        M = sparse(None)
        M.load('./testData/test.haps.sm')
        self.assertEqual(M.row_num(0), 3)
        self.assertEqual(M.row_num(1), 0)
        self.assertEqual(M.row_num(8), 1)

    def test_read_legend(self):
        legend_header, legend = read_legend('./testData/test.legend')
        self.assertEqual(legend_header, ['id','position','a0','a1'])
        self.assertEqual(legend[0], {'id':'19:14492336_A_G', 'position':'1', 'a0':'A', 'a1':'G'})

    def test_assign_bins(self):
        legend_header, legend = read_legend('./testData/test.legend')
        M = sparse(None)
        M.load('./testData/test.haps.sm')
        bins = read_expected('./testData/testBins.txt')

        bin_h = assign_bins(M, bins, legend, False, False, False)
        correct_assignment = {
            0 : [2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15, 18, 21, 22, 23, 24, 25, 26, 29],
            1 : [9, 16, 17, 19, 20, 27, 28, 30],
            2 : [0, 5]
        }
        self.assertEqual(bin_h, correct_assignment)

    def test_verify_legend(self):
        M = sparse(None)
        M.load('./testData/test.haps.sm')
        legend_header, legend = read_legend('./testData/test.nomatch.legend')
        with self.assertRaises(Exception):
            verify_legend(legend, legend_header, M, False)
        with self.assertRaises(Exception):
            verify_legend(legend, legend_header, M, True)

    def test_prune_bins(self):
        M = sparse(None)
        M.load('./testData/test.haps.sm')
        legend_header, legend = read_legend('./testData/test.legend')
        bins = read_expected('./testData/testBins.txt')
        random.seed(123)
        bin_h = bin_h = assign_bins(M, bins, legend, False, False, False)
        prune_bins(bin_h, bins, [], M, 10, 10)
        all_kept_rows = get_all_kept_rows(bin_h, [], False, False, False, False, None)
        true_kept_rows = [0, 4, 5, 8, 9, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30]
        self.assertEqual(all_kept_rows, true_kept_rows)


if __name__ == '__main__':
    unittest.main()