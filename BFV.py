# BFV

from poly import *
from rns import *
from base_converter import *

class BFV:
    # Definitions
    # f(x) = x^n + 1 where n=power-of-two
    # Each HomMult is followed by Relin
    # assumed that we have one q

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
    # -- np (NTT parameters: [w,w_inv,psi,psi_inv])
    # (Generated with parameters)
    # -- sk
    # -- pk
    # -- rlk

    # Extra parameters
    # -- L (max mult depth)
    # -- S (security level)

    def __init__(self, n, q, t, Q, B, qnp, Qnp):
        self.n = n
        self.q = q
        self.Q = Q
        self.t = t
        self.B = B
        self.qnp= qnp # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.Qnp= Qnp # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.rns = RNS(n,q,t, Q,qnp,Qnp)

        ### Calculations ##
        self.q_total = self.rns.base_q.base_prod
        self.q_mod_t = self.q_total % self.t

        self.coeff_div_plain_modulus = self.q_total // self.t
        self.upper_half_increment = self.q_total % self.t

        self.coeff_div_plain_modulus_rns = self.rns.base_q.decompose(self.coeff_div_plain_modulus)
        self.upper_half_increment_rns = self.rns.base_q.decompose(self.upper_half_increment)


    #
    def SecretKeyGen(self):
        """
        sk <- R_2
        """
        s = Poly(self.n,self.q,self.qnp)
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
        a, e = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)
        # a.randomize(self.q)
        a.randomizeV2(self.q)
        e.randomize(self.B)
        pk0 = -(a*self.sk + e)
        pk1 = a
        self.pk = [pk0,pk1]
    #
    def EvalKeyGenV1(self, T):
        self.T = T
        self.l = int(math.floor(math.log(self.q,self.T)))

        rlk1 = []

        sk2 = Poly(self.n,self.q,self.qnp)
        sk2 = (self.sk * self.sk)

        for i in range(self.l+1):
            ai   , ei    = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)
            ai.randomize(self.q)
            ei.randomize(self.B)

            Ts2   = Poly(self.n,self.q,self.qnp)
            Ts2.F = [(self.T**i)*j % self.q for j in sk2.F]

            rlki0 = Ts2 - (ai*self.sk + ei)
            rlki1 = ai

            rlk1.append([rlki0,rlki1])

        self.rlk1 = rlk1
    #
    def EvalKeyGenV2(self, p):
        rlk2 = []
        self.rlk2 = rlk2
        pass
    #
    def Encryption(self, m):
        delta = int(math.floor(self.q/self.t))

        # numerator = (m * self.q_mod_t) + ((self.t + 1) >> 1)
        # fix = numerator // self.t

        u, e1, e2 = Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp), Poly(self.n,self.q,self.qnp)

        u.randomize(2)
        e1.randomize(self.B)
        e2.randomize(self.B)

        c0 = self.pk[0]*u + e1
        c1 = self.pk[1]*u + e2

        # md = Poly(self.n, self.q, self.qnp)
        #
        # md.F = [(delta*x) % self.q for x in m.F]
        temp = [0] * len(self.q)
        for i in range(len(self.q)):
            temp[i] = (self.coeff_div_plain_modulus_rns[i] * m + fix + c0.F[0][i]) % self.q[i]
        c0.F[0] = temp
        # c0 = c0 + md
        return [c0,c1]
    #
    def Decryption(self, ct):
        m = ct[1]*self.sk + ct[0]
        # m.F = [((self.t*x)/self.q) for x in m.F]
        # m = round(m)
        # m = m % self.t
        m = self.rns.decrypt_scale_and_round(m)
        return m
    #
    def DecryptionV2(self, ct):
        sk2 = Poly(self.n,self.q,self.qnp)
        sk2 = (self.sk * self.sk)
        m = (ct[1]*self.sk + ct[0])
        m = (ct[2]*sk2 + m)
        m.F = [((self.t*x)/self.q) for x in m.F]
        m = round(m)
        m = m % self.t
        return m
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
    #
    def __str__(self):
        str = ""
        str = str + "n  : {}\n".format(self.n)
        str = str + "q  : {}\n".format(self.q)
        str = str + "Q  : {}\n".format(self.Q)
        str = str + "t  : {}\n".format(self.t)
        str = str + "B  : {}\n".format(self.B)
        str = str + "qnp: {}\n".format(self.qnp)
        str = str + "Qnp: {}".format(self.Qnp)
        return str
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



        print("ct0[0]:{}".format(ct0[0]))
        print("ct0[1]:{}".format(ct0[1]))
        print("ct1[0]:{}".format(ct1[0]))
        print("ct1[1]:{}".format(ct1[1]))
        # q --> Q
        ct00  = Poly(self.n,self.Q,self.Qnp)
        ct00.F= [x % self.Q for x in ct0[0].F]
        ct01  = Poly(self.n,self.Q,self.Qnp)
        ct01.F= [x % self.Q for x in ct0[1].F]
        ct10  = Poly(self.n,self.Q,self.Qnp)
        ct10.F= [x % self.Q for x in ct1[0].F]
        ct11  = Poly(self.n,self.Q,self.Qnp)
        ct11.F= [x % self.Q for x in ct1[1].F]

        c0 = (ct00 * ct10)
        c1 = (ct00 * ct11) + (ct01 * ct10)
        c2 = (ct01 * ct11)

        print("c0:{}".format(c0))
        print("c1:{}".format(c1))
        print("c2:{}".format(c2))

        c0.F = [((self.t*x)/self.q) for x in c0.F]
        c1.F = [((self.t*x)/self.q) for x in c1.F]
        c2.F = [((self.t*x)/self.q) for x in c2.F]

        c0 = round(c0)
        c1 = round(c1)
        c2 = round(c2)

        # Q --> q

        c0 = (c0 % self.q)
        c1 = (c1 % self.q)
        c2 = (c2 % self.q)

        print("c0:{}".format(c0))
        print("c1:{}".format(c1))
        print("c2:{}".format(c2))

        return [c0,c1,c2]
#
