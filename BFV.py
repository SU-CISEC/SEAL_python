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

    def __init__(self, n, t, B, q, qnp, b, bnp):
        self.n = n
        self.t = t
        self.q = q
        self.b = b
        self.B = B

        # RelinKeys
        self.T = 0
        self.l = 0
        self.p = 0

        # for RNS it changed to [[w,w_inv,psi,psi_inv],[w,w_inv,psi,psi_inv],...]
        self.qnp= qnp # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.bnp = bnp
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
        str = str + "B  : {}\n".format(self.B)
        str = str + "p  : {}\n".format(self.p)
        str = str + "T  : {}\n".format(self.T)
        str = str + "l  : {}\n".format(self.l)

        return str
    #
    def SecretKeyGen(self):
        """
        sk <- R_2
        """
        s = RNSPoly(self.n,self.q,self.qnp)
        s.randomize(2)
        self.sk = s
    #
    def PublicKeyGen(self):
        """
        a <- R_q
        e <- X
        pk[0] <- (-(a*sk)+e) mod q
        pk[1] <- a
        """
        a, e = RNSPoly(self.n,self.q,self.qnp), RNSPoly(self.n,self.q,self.qnp)
        a.randomize(self.q)
        e.randomize(self.B)
        pk0 = -(a*self.sk + e)
        pk1 = a
        self.pk = [pk0,pk1]
    #
    def EvalKeyGenV1(self, T):
        self.T = T
        self.l = int(math.floor(math.log(self.q,self.T)))

        rlk1 = []

        sk2 = (self.sk * self.sk)

        for i in range(self.l+1):
            ai   , ei    = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)
            ai.randomize(self.q)
            ei.randomize(self.B)

            Ts2   = Poly(self.n,self.q,self.qnp)
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
    def Encryption(self, m):
        """
        delta = floor(q/t)

        u  <- random polynomial from R_2
        e1 <- random polynomial from R_B
        e2 <- random polynomial from R_B

        c0 <- pk0*u + e1 + m*delta
        c1 <- pk1*u + e2
        """
        delta = int(math.floor(self.q/self.t))

        u, e1, e2 = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)

        u.randomize(2)
        e1.randomize(self.B)
        e2.randomize(self.B)

        md = Poly(self.n,self.q,self.qnp)
        md.F = [(delta*x) % self.q for x in m.F]

        c0 = self.pk[0]*u + e1
        c0 = c0 + md
        c1 = self.pk[1]*u + e2

        return [c0,c1]
    #
    def Decryption(self, ct):
        """
        ct <- c1*s + c0
        ct <- floot(ct*(t/q))
        m <- [ct]_t
        """
        m = ct[1]*self.sk + ct[0]
        m.F = [((self.t*x)/self.q) for x in m.F]
        m = round(m)
        m = m % self.t
        mr = Poly(self.n,self.t,self.qnp)
        mr.F = m.F
        mr.inNTT = m.inNTT
        return mr
    #
    def DecryptionV2(self, ct):
        """
        ct <- c2*s^2 + c1*s + c0
        ct <- floot(ct*(t/q))
        m <- [ct]_t
        """
        sk2 = (self.sk * self.sk)
        m = ct[0]
        m = (m + (ct[1]*self.sk))
        m = (m + (ct[2]*sk2))
        m.F = [((self.t * x) / self.q) for x in m.F]
        m = round(m)
        m = m % self.t
        mr = Poly(self.n,self.t,self.qnp)
        mr.F = m.F
        mr.inNTT = m.inNTT
        return mr
    #
    def RelinearizationV1(self,ct):
        c0 = ct[0]
        c1 = ct[1]
        c2 = ct[2]

        # divide c2 into base T
        c2i = []

        c2q = Poly(self.n,self.q,self.qnp)
        c2q.F = [x for x in c2.F]

        for i in range(self.l+1):
            c2r = Poly(self.n,self.q,self.qnp)

            for j in range(self.n):
                qt = int(c2q.F[j]/self.T)
                rt = c2q.F[j] - qt*self.T

                c2q.F[j] = qt
                c2r.F[j] = rt

            c2i.append(c2r)

        c0r = Poly(self.n,self.q,self.qnp)
        c1r = Poly(self.n,self.q,self.qnp)
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
    def Encode(self,m): # binary encode
        mr = Poly(self.n,self.t)
        mt = m
        for i in range(self.n):
            mr.F[i] = (mt % self.t)
            mt      = (mt // self.t)
        return mr
    #
    def Decode(self,m): # binary decode
        mr = 0
        for i,c in enumerate(m.F):
            mr = (mr + (c * pow(self.t,i)))# % self.t
        return mr
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
    def HomomorphicMultiplication(self, ct0, ct1):
        ct00 = ct0[0]
        ct01 = ct0[1]
        ct10 = ct1[0]
        ct11 = ct1[1]

        r0 = RefPolMulv2(ct00.F,ct10.F)
        r1 = RefPolMulv2(ct00.F,ct11.F)
        r2 = RefPolMulv2(ct01.F,ct10.F)
        r3 = RefPolMulv2(ct01.F,ct11.F)

        c0 = [x for x in r0]
        c1 = [x+y for x,y in zip(r1,r2)]
        c2 = [x for x in r3]

        c0 = [((self.t * x) / self.q) for x in c0]
        c1 = [((self.t * x) / self.q) for x in c1]
        c2 = [((self.t * x) / self.q) for x in c2]

        c0 = [(round(x) % self.q) for x in c0]
        c1 = [(round(x) % self.q) for x in c1]
        c2 = [(round(x) % self.q) for x in c2]

        # Move to regular modulus
        r0 = Poly(self.n,self.q,self.qnp)
        r1 = Poly(self.n,self.q,self.qnp)
        r2 = Poly(self.n,self.q,self.qnp)

        r0.F = c0
        r1.F = c1
        r2.F = c2

        return [r0,r1,r2]
    # SEAL Encryption function
    def EncryptionSEAL(self, m):
        delta = int(math.floor(self.q/self.t))

        u, e1, e2 = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)

        u.randomize(2)
        e1.randomize(self.B)
        e2.randomize(self.B)

        md = Poly(self.n,self.q,self.qnp)
        md.F = [_ for _ in m.F]

        for i,c in enumerate(md.F):
            if(c >= (self.t >> 1)):
                temp = delta * c

                temp_h = temp >> 64
                temp_l = temp & 0xFFFFFFFFFFFFFFFF

                temp_l = temp_l + 1
                temp_h = temp_h + (temp_l >> 64)

                temp = (temp_h << 64) + temp_l

                temp = temp % self.q

                md.F[i] = temp
            else:
                md.F[i] = (c * delta) % self.q

        c0 = self.pk[0]*u + e1
        c0 = c0 + md
        c1 = self.pk[1]*u + e2

        return [c0,c1]
    # SEAL Decryption function
    def DecryptionSEAL(self, ct, g):
        """
        ct <- c1*s + c0
        st,sg <- fastbconv(ct)
        m <- (st-sg) and scaling
        m <- m* (g^-1 mod t)
        """
        m = ct[1]*self.sk + ct[0]

        g_div_2 = (g>>1)
        gt      = (g*self.t)%self.q
        q_inv_t = modinv((-self.q)%self.t, self.t)
        q_inv_g = modinv((-self.q)%g, g)
        g_inv_t = modinv(g, self.t)

        m.F = [(x*gt) % self.q for x in m.F]

        st, sg = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)

        st.F = [(x*q_inv_t) % self.t for x in m.F]
        sg.F = [(x*q_inv_g) % g      for x in m.F]

        for i in range(self.n):
            m.F[i] = (st.F[i] - sg.F[i]) % self.t

            if(sg.F[i] > g_div_2):
                m.F[i] = (m.F[i] + g) % self.t
            else:
                m.F[i] = m.F[i]

        m.F = [(x*g_inv_t) % self.t for x in m.F]

        mr = Poly(self.n,self.t,self.qnp)
        mr.F = m.F
        mr.inNTT = m.inNTT

        return mr
#

    def EncryptionSEAL_RNS(self, m, m_value):
        # delta = int(math.floor(self.q/self.t))
        # Decompose q_mod_t into rns

        u, e1, e2 = RNSPoly(self.n,self.q,self.qnp), RNSPoly(self.n,self.q,self.qnp), RNSPoly(self.n,self.q,self.qnp)

        u.randomize(2)
        e1.randomize(self.B)
        e2.randomize(self.B)

        numerator = (m_value * self.q_mod_t) + ((self.t + 1) >> 1)
        fix = numerator // self.t

        md = RNSPoly(self.n,self.q[:-1],self.qnp[-1])
        # Copy message m into RNS poly
        # md.copy_poly(m)

        c0 = self.pk[0] * u + e1
        c1 = self.pk[1] * u + e2
        c0 = self.rns.divide_and_round_q_last_inplace(c0)
        c1 = self.rns.divide_and_round_q_last_inplace(c1)

        for i in range(len(self.q) - 1):
            md[i][0] = ((m_value * self.qi_div_t_rns[i]) + fix)

        c0.drop_last_poly()
        c1.drop_last_poly()
        c0 = c0 + md

        return [c0, c1]

    def EncodeSEAL(self,m): # hex encode
        hex_poly = str(m)
        size = len(hex_poly)
        mr = Poly(self.n, self.t)
        value = 0
        for i in range(len(hex_poly)):
            value += int(hex_poly[size - i - 1]) * pow(16, i)
            mr.F[size - i - 1] = int(hex_poly[i])
        return mr, value
    #
    def DecodeSEAL(self,m): # binary decode
        mr = 0
        for i,c in enumerate(m.F):
            mr = (mr + (c * pow(16, i))) % self.t
        return mr

# SEAL Decryption RNS function
    def DecryptionSEAL_RNS(self, ct):
        """
        ct <- c1*s + c0
        st,sg <- fastbconv(ct)
        m <- (st-sg) and scaling
        m <- m* (g^-1 mod t)
        """
        self.sk.drop_last_poly()
        m = ct[1]*self.sk + ct[0]
        print('After dot prod with sk:' + str(m))
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

