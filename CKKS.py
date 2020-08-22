from poly import *


class CKKS:
    def __init__(self, n, t, q, ntt_q):
        self.n = n
        self.t = t
        self.q = q
        self.ntt_q = ntt_q
        self.slots = self.n >> 1


    def encode(self):
        pass

    def decode(self):
        pass

    def SecretKeyGen(self):
        """
        sk <- R_2
        """
        s = RNSPoly(self.n, self.q, self.qnp)
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
        a, e = RNSPoly(self.n, self.q, self.qnp), RNSPoly(self.n, self.q, self.qnp)
        a.randomize(self.q)
        e.randomize(self.B)
        pk0 = a * self.sk
        pk0 = pk0 + e
        pk0 = -pk0
        pk1 = a
        self.pk = [pk0, pk1]

    def encrypt(self):
        pass

    def decrypt(self):
        pass

    def multiply(self):
        pass

    def add(self):
        pass

    def sub(self):
        pass
