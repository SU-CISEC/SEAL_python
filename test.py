import unittest
from rns import *
from poly import *
from ntt import *
from BFV import *
from ntt import *

class BFVTest(unittest.TestCase):
    def test_encrypt_decrypt(self):
        t = 1 << 6
        n = 256

        q = [1099511590913, 1099511592961, 1099511603713]
        k = len(q)
        q_root = [12044354199, 168913082, 694658335]

        b = [2305843009213317121, 2305843009213243393]
        l = len(b)
        b_root = [829315415491244, 32973993658837]

        # Determine B (bound of distribution X)
        sigma = 3.1
        B = int(10 * sigma)

        # Determine T, p (relinearization)
        T = 256
        p = 16

        Evaluator = BFV(n, t, q, b, q_root, b_root)

        # Generate Keys
        Evaluator.SecretKeyGen()
        Evaluator.PublicKeyGen()

        # Encode message
        m_poly, m_len = Evaluator.EncodeInt(0x12345678)
        encrypted = Evaluator.EncryptionSEAL_RNS(m_poly, m_len)

        decrypted = Evaluator.DecryptionSEAL_RNS(encrypted)
        print(decrypted)

        result = Evaluator.DecodeInt(decrypted, m_len)

        self.assertEqual(0x12345678, result)

    def test_bfv_mult(self):
        t = 1 << 6
        n = 256

        q = [1099511590913, 1099511592961, 1099511603713]
        k = len(q)
        q_root = [12044354199, 168913082, 694658335]
        q_ntt = [generate_NTT_tables(n, q[i], q_root[i]) for i in range(k)]

        b = [2305843009213317121, 2305843009213243393]
        l = len(b)
        b_root = [829315415491244, 32973993658837]
        b_ntt = [generate_NTT_tables(n, b[i], b_root[i]) for i in range(l)]

        # Determine B (bound of distribution X)
        sigma = 3.1
        B = int(10 * sigma)

        # Determine T, p (relinearization)
        T = 256
        p = 16

        Evaluator = BFV(n, t, q, b, q_ntt, b_ntt)

        # Generate Keys
        Evaluator.SecretKeyGen()
        Evaluator.PublicKeyGen()

        # Encode message
        m_poly, m_len = Evaluator.EncodeInt(10)
        encrypted = Evaluator.EncryptionSEAL_RNS(m_poly, m_len)

        m_poly2, m_len2 = Evaluator.EncodeInt(5)
        encrypted2 = Evaluator.EncryptionSEAL_RNS(m_poly2, m_len2)

        mul_enc = Evaluator.HomomorphicMultiplication_RNS(encrypted, encrypted2)

        decrypted = Evaluator.DecryptionSEAL_RNS(mul_enc)

        result = Evaluator.DecodeInt(decrypted, m_len + m_len2)

        self.assertEqual(50, result)






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

    def test_decrypt_scale_and_round(self):
        n = 2
        t = 3
        q = [5, 7]
        b = [2305843009213317121, 2305843009213243393]


        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 35
        rns_poly[0][1] = 70
        rns_poly[1][0] = 35
        rns_poly[1][1] = 70

        result = test_rns.decrypt_scale_and_round(rns_poly)
        self.assertEqual(0, result[0])
        self.assertEqual(0, result[1])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 29
        rns_poly[0][1] = 30 + 35
        rns_poly[1][0] = 29
        rns_poly[1][1] = 30 + 35

        result = test_rns.decrypt_scale_and_round(rns_poly)
        self.assertEqual(2, result[0])
        self.assertEqual(0, result[1])

    def test_divide_and_round_q_last_inplace(self):
        """Case 1"""
        n = 2
        t = 3
        q = [13, 7]
        b = [2305843009213317121, 2305843009213243393]

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 1
        rns_poly[0][1] = 2
        rns_poly[1][0] = 1
        rns_poly[1][1] = 2
        result = test_rns.divide_and_round_q_last_inplace(rns_poly)

        self.assertEqual(0, result[0][0])
        self.assertEqual(0, result[0][1])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 12
        rns_poly[0][1] = 11
        rns_poly[1][0] = 4
        rns_poly[1][1] = 3
        result = test_rns.divide_and_round_q_last_inplace(rns_poly)

        self.assertEqual(4, result[0][0])
        self.assertEqual(3, result[0][1])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 6
        rns_poly[0][1] = 2
        rns_poly[1][0] = 5
        rns_poly[1][1] = 1
        result = test_rns.divide_and_round_q_last_inplace(rns_poly)

        self.assertEqual(3, result[0][0])
        self.assertEqual(2, result[0][1])

        """Case 2"""
        n = 2
        t = 1
        q = [3, 5, 7, 11]
        b = [2305843009213317121, 2305843009213243393]

        test_rns = RNS(n, t, q, b, qnp=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], bnp=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 1
        rns_poly[0][1] = 2
        rns_poly[1][0] = 1
        rns_poly[1][1] = 2
        rns_poly[2][0] = 1
        rns_poly[2][1] = 2
        rns_poly[3][0] = 1
        rns_poly[3][1] = 2
        result = test_rns.divide_and_round_q_last_inplace(rns_poly)

        self.assertEqual(0, result[0][0])
        self.assertEqual(0, result[0][1])
        self.assertEqual(0, result[1][0])
        self.assertEqual(0, result[1][1])
        self.assertEqual(0, result[2][0])
        self.assertEqual(0, result[2][1])

        # Next a case with non-trivial rounding; array is (60, 70)
        rns_poly = RNSPoly(n, q)
        rns_poly[0][0] = 0
        rns_poly[0][1] = 1
        rns_poly[1][0] = 0
        rns_poly[1][1] = 0
        rns_poly[2][0] = 4
        rns_poly[2][1] = 0
        rns_poly[3][0] = 5
        rns_poly[3][1] = 4
        result = test_rns.divide_and_round_q_last_inplace(rns_poly)

        # We get only approximate result in this case
        self.assertTrue((3 + 2 - result[0][0]) % 3 <= 1)
        self.assertTrue((3 + 0 - result[0][1]) % 3 <= 1)
        self.assertTrue((5 + 0 - result[1][0]) % 5 <= 1)
        self.assertTrue((5 + 1 - result[1][1]) % 5 <= 1)
        self.assertTrue((7 + 5 - result[2][0]) % 7 <= 1)
        self.assertTrue((7 + 6 - result[2][1]) % 7 <= 1)



