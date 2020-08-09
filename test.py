import unittest
from rns import *
from poly import *


class RNSTest(unittest.TestCase):

    def test_fastbconv_m_tilde(self):
        n = 2
        t = 1
        q = [3, 5]
        b = [2305843009213317121, 2305843009213243393]
        m_tilde = 1 << 32

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])
        rns_poly = RNSPoly(n, q)

        rns_poly[0][0] = 1
        rns_poly[0][1] = 1
        rns_poly[1][0] = 2
        rns_poly[1][1] = 2
        result = test_rns.fastbconv_m_tilde(rns_poly)

        temp = ((2 * m_tilde) % 3) * 5 + ((4 * m_tilde) % 5) * 3;

        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[0], result[0][0])
        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[0], result[0][1])
        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[1], result[1][0])
        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[1], result[1][1])
        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[2], result[2][0])
        self.assertEqual(temp % test_rns.base_Bsk_m_tilde[2], result[2][1])

    def test_sm_mrq(self):
        n = 2
        t = 1
        q = [3, 5]
        b = [2305843009213317121, 2305843009213243393]
        m_tilde = 1 << 32

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly = RNSPoly(n, test_rns.base_Bsk_m_tilde.rns_base)
        rns_poly[0][0] = m_tilde
        rns_poly[0][1] = 2 * m_tilde
        rns_poly[1][0] = m_tilde
        rns_poly[1][1] = 2 * m_tilde
        rns_poly[2][0] = m_tilde
        rns_poly[2][1] = 2 * m_tilde
        rns_poly[3][0] = 0
        rns_poly[3][1] = 0

        # This should simply get rid of the m_tilde factor
        result = test_rns.sm_mrq(rns_poly)

        self.assertEqual(1, result[0][0])
        self.assertEqual(2, result[0][1])
        self.assertEqual(1, result[1][0])
        self.assertEqual(2, result[1][1])
        self.assertEqual(1, result[2][0])
        self.assertEqual(2, result[2][1])

        rns_poly = RNSPoly(n, test_rns.base_Bsk_m_tilde.rns_base)
        rns_poly[0][0] = 15
        rns_poly[0][1] = 30
        rns_poly[1][0] = 15
        rns_poly[1][1] = 30
        rns_poly[2][0] = 15
        rns_poly[2][1] = 30
        rns_poly[3][0] = 15
        rns_poly[3][1] = 30

        # Next add a multiple of q to the input and see if it is reduced properly
        result = test_rns.sm_mrq(rns_poly)

        self.assertEqual(0, result[0][0])
        self.assertEqual(0, result[0][1])
        self.assertEqual(0, result[1][0])
        self.assertEqual(0, result[1][1])
        self.assertEqual(0, result[2][0])
        self.assertEqual(0, result[2][1])

        # Now with a multiple of m_tilde + multiple of q

        rns_poly = RNSPoly(n, test_rns.base_Bsk_m_tilde.rns_base)
        rns_poly[0][0] = m_tilde + 15
        rns_poly[0][1] = m_tilde + 30
        rns_poly[1][0] = 2 * m_tilde + 15
        rns_poly[1][1] = m_tilde + 30
        rns_poly[2][0] = 2 * m_tilde + 15
        rns_poly[2][1] = 2 * m_tilde + 30
        rns_poly[3][0] = 2 * m_tilde + 15
        rns_poly[3][1] = 2 * m_tilde + 30

        result = test_rns.sm_mrq(rns_poly)

        self.assertEqual(1, result[0][0])
        self.assertEqual(1, result[0][1])
        self.assertEqual(2, result[1][0])
        self.assertEqual(1, result[1][1])
        self.assertEqual(2, result[2][0])
        self.assertEqual(2, result[2][1])

    def test_fast_floor(self):
        n = 2
        t = 1
        q = [3, 5]
        b = [2305843009213317121, 2305843009213243393]
        m_tilde = 1 << 32

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly_q = RNSPoly(n,q)
        rns_poly_bsk = RNSPoly(n,test_rns.base_Bsk.rns_base)
        rns_poly_q[0][0] = 15
        rns_poly_q[0][1] = 30
        rns_poly_q[1][0] = 15
        rns_poly_q[1][1] = 30

        rns_poly_bsk[0][0] = 15
        rns_poly_bsk[0][1] = 30
        rns_poly_bsk[1][0] = 15
        rns_poly_bsk[1][1] = 30
        rns_poly_bsk[2][0] = 15
        rns_poly_bsk[2][1] = 30

        result = test_rns.fast_floor(rns_poly_q, rns_poly_bsk)

        self.assertEqual(1, result[0][0])
        self.assertEqual(2, result[0][1])
        self.assertEqual(1, result[1][0])
        self.assertEqual(2, result[1][1])
        self.assertEqual(1, result[2][0])
        self.assertEqual(2, result[2][1])

        rns_poly_q = RNSPoly(n, q)
        rns_poly_bsk = RNSPoly(n, test_rns.base_Bsk.rns_base)
        rns_poly_q[0][0] = 21
        rns_poly_q[0][1] = 32
        rns_poly_q[1][0] = 21
        rns_poly_q[1][1] = 32

        rns_poly_bsk[0][0] = 21
        rns_poly_bsk[0][1] = 32
        rns_poly_bsk[1][0] = 21
        rns_poly_bsk[1][1] = 32
        rns_poly_bsk[2][0] = 21
        rns_poly_bsk[2][1] = 32

        # The result is not exact but differs at most by 1
        result = test_rns.fast_floor(rns_poly_q, rns_poly_bsk)

        self.assertTrue((1 - result[0][0]) <= 1)
        self.assertTrue((2 - result[0][1]) <= 1)
        self.assertTrue((1 - result[1][0]) <= 1)
        self.assertTrue((2 - result[1][1]) <= 1)
        self.assertTrue((1 - result[2][0]) <= 1)
        self.assertTrue((2 - result[2][1]) <= 1)

    def test_fastbconv_sk(self):
        n = 2
        t = 1
        q = [3, 5]
        b = [2305843009213317121, 2305843009213243393]
        m_tilde = 1 << 32

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly_bsk = RNSPoly(n, test_rns.base_Bsk.rns_base)

        rns_poly_bsk[0][0] = 1
        rns_poly_bsk[0][1] = 2
        rns_poly_bsk[1][0] = 1
        rns_poly_bsk[1][1] = 2
        rns_poly_bsk[2][0] = 1
        rns_poly_bsk[2][1] = 2

        result = test_rns.fastbconv_sk(rns_poly_bsk)

        self.assertEqual(1, result[0][0])
        self.assertEqual(2, result[0][1])
        self.assertEqual(1, result[1][0])
        self.assertEqual(2, result[1][1])




    def test_fast_convert_array(self):

        """############CASE 1############"""
        in_base = RNSBase([2, 3])
        out_base = RNSBase([2])

        bct = BaseConverter(in_base, out_base)

        rns_poly = RNSPoly(3, [2, 3])

        rns_poly[0][0] = 0
        rns_poly[0][1] = 1
        rns_poly[0][2] = 0
        rns_poly[1][0] = 0
        rns_poly[1][1] = 1
        rns_poly[1][2] = 2

        result = bct.fast_convert_array(rns_poly)

        self.assertEqual(0, result[0][0])
        self.assertEqual(1, result[0][1])
        self.assertEqual(0, result[0][2])

        """############CASE 2############"""
        in_base = RNSBase([2, 3])
        out_base = RNSBase([2, 3])

        bct = BaseConverter(in_base, out_base)

        rns_poly = RNSPoly(3, [2, 3])

        rns_poly[0][0] = 1
        rns_poly[0][1] = 1
        rns_poly[0][2] = 0
        rns_poly[1][0] = 1
        rns_poly[1][1] = 2
        rns_poly[1][2] = 2

        result = bct.fast_convert_array(rns_poly)

        self.assertEqual(1, result[0][0])
        self.assertEqual(1, result[0][1])
        self.assertEqual(0, result[0][2])
        self.assertEqual(1, result[1][0])
        self.assertEqual(2, result[1][1])
        self.assertEqual(2, result[1][2])

        """############CASE 3############"""
        in_base = RNSBase([2, 3])
        out_base = RNSBase([3, 4, 5])

        bct = BaseConverter(in_base, out_base)

        rns_poly = RNSPoly(3, [2, 3])

        rns_poly[0][0] = 0
        rns_poly[0][1] = 1
        rns_poly[0][2] = 1
        rns_poly[1][0] = 0
        rns_poly[1][1] = 1
        rns_poly[1][2] = 2

        result = bct.fast_convert_array(rns_poly)

        self.assertEqual(0, result[0][0])
        self.assertEqual(1, result[0][1])
        self.assertEqual(2, result[0][2])
        self.assertEqual(0, result[1][0])
        self.assertEqual(3, result[1][1])
        self.assertEqual(1, result[1][2])
        self.assertEqual(0, result[2][0])
        self.assertEqual(2, result[2][1])
        self.assertEqual(0, result[2][2])
