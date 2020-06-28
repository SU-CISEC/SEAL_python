# Poly

from random import randint,sample
from ntt import *


class Poly:
    def __init__(self, n, q, np):
        self.n = n
        self.q = q  # q = [q0,q1,....]
        self.np = np  # NTT parameters: [[w,w_inv,psi,psi_inv],[w,w_inv,psi,psi_inv], .....]
        self.F = [[0] * len(q)] * n
        self.inNTT = False

    #
    def randomize(self, B, domain=False):
        for i in range(len(self.F)):
            self.F[i] = [randint(0, B - 1) for j in range(len(self.q))]
        self.inNTT = domain

    def randomizeV2(self, B, domain=False):
        for i in range(self.n):
            tmp = [0] * len(self.q)
            for j in range(len(B)):
                tmp[j] = randint(0, B[j] - 1)
            self.F[i] = tmp
        self.inNTT = domain

    # # Check reduce of t
    # def reduce(self, t):
    #     self.F = [x % t for x in self.F]
    #     self.q = t
    # #
    def __str__(self):
        tmp = min(self.n, 8)
        str = ""
        # str = str + "n      : {}\n".format(self.n)
        # str = str + "q      : {}\n".format(self.q)
        # str = str + "np     : {}\n".format(self.np)
        # str = str + "In NTT?: {}\n".format(self.inNTT)
        str = str + "First {} coefficients of polynomial: {}".format(tmp, self.F[0:tmp])
        return str

    #
    def __add__(self, b):
        c = Poly(self.n, self.q, self.np)
        for i in range(self.n):
            for j in range(len(self.q)):
                c.F[i][j] = self.F[i][j] + b.F[i][j] % self.q[j]
        c.inNTT = self.inNTT
        return c

    #
    def __sub__(self, b):
        c = Poly(self.n, self.q, self.np)
        for i in range(self.n):
            for j in range(len(self.q)):
                c.F[i][j] = self.F[i][j] - b.F[i][j] % self.q[j]
        c.inNTT = self.inNTT
        return c

    #
    def __mul__(self, b):
        """
        Assuming both inputs in POL/NTT domain
        """
        c = Poly(self.n, self.q, self.np)
        if self.inNTT == True and b.inNTT == True:
            for i in range(self.n):
                for j in range(len(self.q)):
                    c.F[i][j] = self.F[i][j] * b.F[i][j] % self.q[j]
            c.inNTT = True
        else:
            for j in range(len(self.q)):
                # x1=self*psi, x2=b*psi
                # x1n = NTT(x1,w), x2n = NTT(x2,w)
                # x3n = x1n*x2n
                # x3 = INTT(x3n,w_inv)
                # c = x3*psi_inv
                psi, psi_inv, W, W_inv = self.np[j]

                row_a = [row[j] for row in self.F]
                s_p = [(x * (psi ** pwr)) % self.q[j] for pwr, x in enumerate(row_a)]
                row_b = [row[j] for row in b.F]
                b_p = [(x * (psi ** pwr)) % self.q[j] for pwr, x in enumerate(row_b)]
                s_n = NTT(s_p, W, self.q[j])
                b_n = NTT(b_p, W, self.q[j])
                sb_n = [(x * y) % self.q[j] for x, y in zip(s_n, b_n)]
                sb_p = INTT(sb_n, W_inv, self.q[j])

                sb = [(x * (psi_inv ** pwr)) % self.q[j] for pwr, x in enumerate(sb_p)]

                for i in range(len(sb)):
                    c.F[i][j] = int(sb[i])
                c.inNTT = True
        return c

    def __getitem__(self, index):
        return self.F[index]

    def __setitem__(self, key, value):
        print(key)
        self.F[key] = value
    #
    # def multiply_by_scalar_coeff(self, scalar, modulus):
    #

    #
    # # !!!!!!!!!!check!!!!!!!!!!!!!!!!
    # def __mod__(self,base):
    #     b = Poly(self.n, base, self.np)
    #
    #     for i in range(self.n):
    #         for j in range(len(self.q)):
    #             b.F[i][j] = self.F[i][j] % self.q[j]
    #     b.F = [x%base for x in self.F]
    #     b.inNTT = self.inNTT
    #     return b
    #
    def __round__(self):
        b = Poly(self.n, self.q, self.np)
        for i in range(len(self.F)):
            for j in range(len(self.q)):
                b.F[i][j] = round(self.F[i][j])
        b.inNTT = self.inNTT
        return b

    #
    def __eq__(self, b):
        if self.n != b.n:
            return False
        else:
            for i in range(self.n):
                for j in range(len(self.q)):
                    if self.F[i][j] != b.F[i][j]:
                        return False
            return True

    #
    def __neg__(self):
        b = Poly(self.n, self.q, self.np)
        for i in range(self.n):
            for j in range(len(self.q)):
                b.F[i][j] = (-self.F[i][j])
        b.inNTT = self.inNTT
        return b

    #
    def toNTT(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            for j in range(len(self.q)):
                row_j = [row[j] for row in self.F]
                row_j = NTT(row_j, self.np[j][0], self.q)
                for i in range(self.n):
                    b.F[i][j] = row_j[i]
            b.inNTT = True
        else:
            b.F = self.F.copy()
            b.inNTT = True
        return b

    #
    def toPOL(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            b.F = self.F.copy()
            b.inNTT = False
        else:
            for j in range(len(self.q)):
                row_j = [row[j] for row in self.F]
                row_j = INTT(row_j, self.np[j][1], self.q)
                for i in range(self.n):
                    b.F[i][j] = row_j[i]
            b.inNTT = False
        return b
#
