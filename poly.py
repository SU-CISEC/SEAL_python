# Poly

from random import randint
from ntt import *


class RNSPoly:
    def __init__(self,n,RNS_mod,ntt_tables=None):
        self.n = n
        self.k = len(RNS_mod)
        if ntt_tables is not None:
            self.rns_poly = [Poly(n,RNS_mod[i],ntt_tables[i]) for i in range(self.k)]
        else:
            self.rns_poly = [Poly(n, RNS_mod[i]) for i in range(self.k)]
        self.RNS_mod = RNS_mod
        self.ntt_tables = ntt_tables


    def randomize(self,limit):
        if type(limit) is list:
            # Random numbers for different q
            for i in range(self.k):
                self[i].randomize(limit[i])
        else:
            for i in range(self.k):
                self[i].randomize(limit)

    def append(self, b):
        for i in range(b.k):
            self.rns_poly.append(b[i])
            self.RNS_mod.append(b.RNS_mod[i])
            self.ntt_tables.append(b.ntt_tables[i])
        self.k += b.k

    def get_poly(self,key):
        poly = Poly(self.n,self.RNS_mod[key])
        poly.F = self.rns_poly[key].F.copy()
        poly.inNTT = self.rns_poly[key].inNTT
        poly.np = self.ntt_tables[key] if self.ntt_tables is not None else [0,0,0,0]
        return poly

    def __str__(self):
        rns_str = '\n'
        for i in range(self.k):
            rns_str += '[' + str(i) + ']: ' + str(self[i]) + '\n'
        return rns_str

    def __mul__(self, b):
        c = RNSPoly(self.n, self.RNS_mod, self.ntt_tables)
        if type(b) is int:
            for i in range(self.k):
                c[i] = self[i] * b
        else:
            for i in range(self.k):
                c[i] = self[i] * b[i]

        return c

    def __add__(self, b):
        c = RNSPoly(self.n, self.RNS_mod, self.ntt_tables)
        if type(b) is int:
            for i in range(self.k):
                c[i] = self[i] + b
        else:
            for i in range(self.k):
                c[i] = self[i] + b[i]

        return c

    def __sub__(self, b):
        c = RNSPoly(self.n, self.RNS_mod, self.ntt_tables)
        for i in range(self.k):
            c[i] = self[i] - b[i]
        return c

    def __neg__(self):
        c = RNSPoly(self.n, self.RNS_mod, self.ntt_tables)
        for i in range(self.k):
            c[i] = -self[i]
        return c

    def __getitem__(self, key):
        return self.rns_poly[key]

    def __setitem__(self, key, value):
        self.rns_poly[key] = value

    def copy_poly(self, poly):
        c = RNSPoly(self.n, self.RNS_mod, self.ntt_tables)
        for i in range(self.k):
            c[i] = poly % self.RNS_mod[i]

    def drop_last_poly(self):
        self.rns_poly = self.rns_poly[:-1]
        self.ntt_tables = self.ntt_tables[:-1] if self.ntt_tables is not None else None
        self.RNS_mod = self.RNS_mod[:-1]
        self.k = len(self.RNS_mod)



class Poly:
    def __init__(self, n, q, np=[0,0,0,0]):
        self.n = n
        self.q = q
        self.np= np # NTT parameters: [w,w_inv,psi,psi_inv]
        self.F = [0]*n
        self.inNTT = False
    #
    def randomize(self, B, domain=False):
        if B == 2:
            # self.F = [randint(-1, B - 1) % self.q for i in range(self.n)]
            self.F = [1 for i in range(self.n)]
            self.inNTT = domain
        else:
            # self.F = [randint(0, B-1) for i in range(self.n)]
            self.F = [1 for i in range(self.n)]
            self.inNTT = domain
    #
    def __str__(self):
        pstr = str(self.F[0])
        tmp = min(self.n,8)

        for i in range(1,tmp):
            pstr = pstr+" + "+str(self.F[i])+"*x^"+str(i)

        if self.n > 8:
            pstr = pstr + " + ..."
        return pstr
    #
    def __add__(self, b):
        if type(b) == int:
            c = Poly(self.n, self.q, self.np)
            c.F = [((x + b) % self.q) for x in self.F]
            return c
        else:
            if self.inNTT != b.inNTT:
                raise Exception("Polynomial Addiditon: Inputs must be in the same domain.")
            elif self.q != b.q:
                raise Exception("Polynomial Addiditon: Inputs must have the same modulus")
            else:
                c = Poly(self.n, self.q, self.np)
                c.F = [(x+y)%self.q for x,y in zip(self.F,b.F)]
                c.inNTT = self.inNTT
                return c
    #
    def __sub__(self, b):
        if type(b) == int:
            c = Poly(self.n, self.q, self.np)
            c.F = [((x - b) % self.q) for x in self.F]
            return c
        else:
            if self.inNTT != b.inNTT:
                raise Exception("Polynomial Subtraction: Inputs must be in the same domain.")
            elif self.q != b.q:
                raise Exception("Polynomial Subtraction: Inputs must have the same modulus")
            else:
                c = Poly(self.n, self.q, self.np)
                c.F = [(x-y)%self.q for x,y in zip(self.F,b.F)]
                c.inNTT = self.inNTT
                return c
    #
    def __mul__(self, b):
        if type(b) == int:
            c = Poly(self.n, self.q, self.np)
            c.F = [((x * b) % self.q) for x in self.F]
            return c
        else:
            if self.inNTT != b.inNTT:
                raise Exception("Polynomial Multiplication: Inputs must be in the same domain.")
            elif self.q != b.q:
                raise Exception("Polynomial Multiplication: Inputs must have the same modulus")
            else:
                """
                Assuming both inputs in POL/NTT domain
                If in NTT domain --> Coeff-wise multiplication
                If in POL domain --> Full polynomial multiplication
                """
                c = Poly(self.n, self.q, self.np)
                if self.inNTT == True and b.inNTT == True:
                    c.F = [((x*y)%self.q) for x,y in zip(self.F,b.F)]
                    c.inNTT = True
                else:
                    # x1=self*psi, x2=b*psi
                    # x1n = NTT(x1,w), x2n = NTT(x2,w)
                    # x3n = x1n*x2n
                    # x3 = INTT(x3n,w_inv)
                    # c = x3*psi_inv

                    w_table    = self.np[0]
                    wv_table   = self.np[1]
                    psi_table  = self.np[2]
                    psiv_table = self.np[3]

                    s_p = [(x*psi_table[pwr])%self.q for pwr,x in enumerate(self.F)]
                    b_p = [(x*psi_table[pwr])%self.q for pwr,x in enumerate(b.F)]
                    s_n = NTT(s_p,w_table,self.q)
                    b_n = NTT(b_p,w_table,self.q)
                    sb_n= [(x*y)%self.q for x,y in zip(s_n,b_n)]
                    sb_p= INTT(sb_n,wv_table,self.q)
                    sb  = [(x*psiv_table[pwr])%self.q for pwr,x in enumerate(sb_p)]

                    c.F = sb
                    c.inNTT = False
                return c
    #
    def __mod__(self,base):
        b = Poly(self.n, base, self.np)
        b.F = [(x%base) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def __round__(self):
        b = Poly(self.n, self.q, self.np)
        b.F = [round(x) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def __eq__(self, b):
        if self.n != b.n:
            return False
        elif self.q != b.q:
            return False
        else:
            for i,j in zip(self.F,b.F):
                if i != j:
                    return False
            return True
    #
    def __neg__(self):
        b = Poly(self.n, self.q, self.np)
        b.F = [((-x) % self.q) for x in self.F]
        b.inNTT = self.inNTT
        return b

    def __getitem__(self, key):
        return self.F[key]

    def __setitem__(self, key, value):
        self.F[key] = value

    def __floordiv__(self, b):
        b = Poly(self.n, self.q, self.np)
        b.F = [(x // b) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def toNTT(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            b.F = NTT(self.F,self.np[0],self.q)
            b.inNTT = True
        else:
            b.F = [x for x in self.F]
            b.inNTT = True
        return b
    #
    def toPOL(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            b.F = [x for x in self.F]
            b.inNTT = False
        else:
            b.F = INTT(self.F,self.np[1],self.q)
            b.inNTT = False
        return b

    def to_rns_poly(self):
        rns_poly = RNSPoly(self.n, [self.q], [self.np])
        rns_poly.k = 1
        rns_poly.rns_poly[0].F = self.F.copy()
        rns_poly.rns_poly[0].inNTT = self.inNTT
        return rns_poly


#
