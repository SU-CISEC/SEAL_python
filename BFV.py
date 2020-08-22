# BFV

from poly import *
from rns import *

class BFV:
    # Definitions
    # Z_q[x]/f(x) = x^n + 1 where n=power-of-two

    # Operations
    # -- SecretKeyGen
    # -- PublicKeyGen
    # -- Encryption
    # -- Decryption
    # -- EvaluationKeyGen
    # -- HomAdd
    # -- HomMult
    # -- RelinV1
    # -- RelinV2

    # Parameters
    # (From outside)
    # -- n (ring size)
    # -- q (ciphertext modulus)
    # -- t (plaintext modulus)
    # -- B (distribution X bound)
    # -- qnp (NTT parameters: [w,w_inv,psi,psi_inv])
    # (Generated with parameters)
    # -- sk
    # -- pk
    # -- rlk

    # Extra parameters
    # -- L (max mult depth)
    # -- S (security level)

    def __init__(self, n, t, q, b, qnp, bnp):
        self.n = n
        self.t = t
        self.q = q
        self.b = b

        # RelinKeys
        self.T = 0
        self.l = 0
        self.p = 0

        # for RNS it changed to [[psi_0],[psi_1],...]
        self.q_ntt_table = qnp # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.b_ntt_table = bnp
        #
        self.sk = []
        self.pk = []
        self.rlk1 = []
        self.rlk2 = []

        # RNS init
        self.rns = RNS(n, t, q, b, qnp, bnp)
        self.rns2 = RNS(n, t, q[:-1], b[:-1], qnp[:-1], bnp[:-1])

        # Calculations #
        self.q_total = self.rns.base_q.base_prod
        self.q_mod_t = self.q_total % self.t
        self.q_div_t = self.q_total // self.t

        self.qi_div_t_rns = self.rns.base_q.decompose(self.q_div_t)
        self.qi_mod_t_rns = self.rns.base_q.decompose(self.q_mod_t)
    #
    def __str__(self):
        str = "\n--- Parameters:\n"
        str = str + "n  : {}\n".format(self.n)
        str = str + "t  : {}\n".format(self.t)
        str = str + "q  : {}\n".format(self.q)
        str = str + "b  : {}\n".format(self.b)
        str = str + "p  : {}\n".format(self.p)
        str = str + "T  : {}\n".format(self.T)
        str = str + "l  : {}\n".format(self.l)

        return str
    #
    def SecretKeyGen(self):
        """
        sk <- R_2
        """
        s = RNSPoly(self.n, self.q, self.q_ntt_table)
        s.randomize(0)
        self.sk = s
    #
    def PublicKeyGen(self):
        """
        a <- R_q
        e <- X
        pk[0] <- (-(a*sk)+e) mod q
        pk[1] <- a
        """
        a, e = RNSPoly(self.n, self.q, self.q_ntt_table), RNSPoly(self.n, self.q, self.q_ntt_table)
        a.randomize(2, True)
        e.randomize(1)
        pk0 = a*self.sk
        pk0 = pk0 + e
        pk0 = -pk0
        pk1 = a
        self.pk = [pk0,pk1]
    #
    def EvalKeyGenV1(self, T):
        self.T = T
        self.l = int(math.floor(math.log(self.q,self.T)))

        rlk1 = []

        sk2 = (self.sk * self.sk)

        for i in range(self.l+1):
            ai   , ei    = Poly(self.n, self.q, self.q_ntt_table), Poly(self.n, self.q, self.q_ntt_table)
            ai.randomize(self.q)
            ei.randomize(self.B)

            Ts2   = Poly(self.n, self.q, self.q_ntt_table)
            Ts2.F = [((self.T**i)*j) % self.q for j in sk2.F]

            rlki0 = Ts2 - (ai*self.sk + ei)
            rlki1 = ai

            rlk1.append([rlki0,rlki1])

        self.rlk1 = rlk1
    #
    def EvalKeyGenV2(self, p):
        self.p = p
        pass
    #
    #
    def RelinearizationV1(self,ct):
        c0 = ct[0]
        c1 = ct[1]
        c2 = ct[2]

        # divide c2 into base T
        c2i = []

        c2q = Poly(self.n, self.q, self.q_ntt_table)
        c2q.F = [x for x in c2.F]

        for i in range(self.l+1):
            c2r = Poly(self.n, self.q, self.q_ntt_table)

            for j in range(self.n):
                qt = int(c2q.F[j]/self.T)
                rt = c2q.F[j] - qt*self.T

                c2q.F[j] = qt
                c2r.F[j] = rt

            c2i.append(c2r)

        c0r = Poly(self.n, self.q, self.q_ntt_table)
        c1r = Poly(self.n, self.q, self.q_ntt_table)
        c0r.F = [x for x in c0.F]
        c1r.F = [x for x in c1.F]

        for i in range(self.l+1):
            c0r = c0r + (self.rlk1[i][0] * c2i[i])
            c1r = c1r + (self.rlk1[i][1] * c2i[i])

        return [c0r,c1r]
    #
    def RelinearizationV2(self,ct):
        pass
    #
    def HomomorphicAddition(self, ct0, ct1):
        ct0_b = ct0[0] + ct1[0]
        ct1_b = ct0[1] + ct1[1]
        return [ct0_b,ct1_b]
    #
    def HomomorphicSubtraction(self, ct0, ct1):
        ct0_b = ct0[0] - ct1[0]
        ct1_b = ct0[1] - ct1[1]
        return [ct0_b,ct1_b]
    #

    def EncryptionSEAL_RNS(self, m_poly, m_len):
        # delta = int(math.floor(self.q/self.t))
        # Decompose q_mod_t into rns

        u, e1, e2 = RNSPoly(self.n, self.q, self.q_ntt_table), RNSPoly(self.n, self.q, self.q_ntt_table), RNSPoly(self.n, self.q, self.q_ntt_table)

        u.randomize(0)
        e1.randomize(1)
        e2.randomize(1)

        c0 = self.pk[0] * u
        c1 = self.pk[1] * u
        c0 = c0 + e1
        c1 = c1 + e2
        c0 = self.rns.divide_and_round_q_last_inplace(c0)
        c1 = self.rns.divide_and_round_q_last_inplace(c1)

        for j in range(m_len):
            numerator = (m_poly[j] * self.q_mod_t) + ((self.t + 1) >> 1)
            fix = numerator // self.t

            for i in range(len(self.q) - 1):
                c0[i][j] = (c0[i][j] + ((m_poly[j] * self.qi_div_t_rns[i]) + fix)) % self.q[i]

        c0.drop_last_poly()
        c1.drop_last_poly()

        return [c0, c1]

    def EncodeInt(self,value):
        m_poly = Poly(self.n, self.t)
        i = 0
        while (value != 0):
            if (value & 1) != 0:
                m_poly[i] = 1
            value >>= 1
            i += 1
        return m_poly, i

    def DecodeInt(self, poly, poly_len):
        coeff_neg_threshold = (self.t + 1) >> 1
        result = 0
        result_is_negative = False
        for i in range(poly_len-1, -1, -1):
            coeff = poly[i]
            result = result << 1
            coeff_is_negative = coeff >= coeff_neg_threshold
            pos_value = coeff

            if coeff_is_negative:
                pos_value = self.t - pos_value

            # Add or subtract-in coefficient.
            if result_is_negative == coeff_is_negative:
                result += pos_value
            else:
                result -= pos_value
                if result < 0:
                    result = -result
                    result_is_negative = not result_is_negative

        return result



# SEAL Decryption RNS function
    def DecryptionSEAL_RNS(self, ct):
        """
        ct <- c1*s + c0
        st,sg <- fastbconv(ct)
        m <- (st-sg) and scaling
        m <- m* (g^-1 mod t)
        """
        self.sk.drop_last_poly()
        mask = ct[1]*self.sk
        m = mask + ct[0]
        mr = self.rns2.decrypt_scale_and_round(m)
        return mr

    def HomomorphicMultiplication_RNS(self, ct0, ct1):
        # Extract RNS_poly
        ct00 = ct0[0]
        ct01 = ct0[1]
        ct10 = ct1[0]
        ct11 = ct1[1]

        # 1 - Do base extension from base q to base (bsk U m_tilde)
        ct00_bsk_mtilde = self.rns2.fastbconv_m_tilde(ct00)
        ct01_bsk_mtilde = self.rns2.fastbconv_m_tilde(ct01)
        ct10_bsk_mtilde = self.rns2.fastbconv_m_tilde(ct10)
        ct11_bsk_mtilde = self.rns2.fastbconv_m_tilde(ct11)

        # 2 - Do small montgomery reduction to get rid of m_tilde
        # It will convert from base (bsk U m_tilde) to base (bsk)
        ct00_bsk = self.rns2.sm_mrq(ct00_bsk_mtilde)
        ct01_bsk = self.rns2.sm_mrq(ct01_bsk_mtilde)
        ct10_bsk = self.rns2.sm_mrq(ct10_bsk_mtilde)
        ct11_bsk = self.rns2.sm_mrq(ct11_bsk_mtilde)

        # 3 - Do dyadic product both base q and base bsk

        # 3(a) - Do for base q
        cq0 = ct00 * ct10
        cq1 = (ct00 * ct11) + (ct01 * ct10)
        cq2 = ct01 * ct11

        # 3(b) - Do for base bsk
        cbsk0 = ct00_bsk * ct10_bsk
        cbsk1 = (ct00_bsk * ct11_bsk) + (ct01_bsk * ct10_bsk)
        cbsk2 = ct01_bsk * ct11_bsk

        # 4 - Do multiply base q and base bsk components by t (plain_modulus)
        cq0_t = cq0 * self.t
        cq1_t = cq1 * self.t
        cq2_t = cq2 * self.t

        cbsk0_t = cbsk0 * self.t
        cbsk1_t = cbsk1 * self.t
        cbsk2_t = cbsk2 * self.t

        # 5 - Do divide by q and floor, producing a result in base Bsk
        cbsk0_q = self.rns2.fast_floor(cq0_t, cbsk0_t)
        cbsk1_q = self.rns2.fast_floor(cq1_t, cbsk1_t)
        cbsk2_q = self.rns2.fast_floor(cq2_t, cbsk2_t)

        # 6 - use Shenoy-Kumaresan method to convert the result to base q and write to encrypted1

        c0 = self.rns2.fastbconv_sk(cbsk0_q)
        c1 = self.rns2.fastbconv_sk(cbsk1_q)
        c2 = self.rns2.fastbconv_sk(cbsk2_q)

        return [c0, c1, c2]

