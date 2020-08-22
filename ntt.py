import math
from helper import *

# --- NTT functions ---

def generate_root_powers(n,q,root):
    psi = root
    psi_table = [1] * n
    for i in range(1, n):
        psi_table[i] = ((psi_table[i - 1] * psi) % q)

    return psi_table
def generate_NTT_tables(n, q, root):
    psi = root
    psiv = modinv(psi, q)
    w = pow(psi, 2, q)
    wv = modinv(w, q)

    w_table = [1] * n
    wv_table = [1] * n
    psi_table = [1] * n
    psiv_table = [1] * n
    for i in range(1, n):
        w_table[i] = ((w_table[i - 1] * w) % q)
        wv_table[i] = ((wv_table[i - 1] * wv) % q)
        psi_table[i] = ((psi_table[i - 1] * psi) % q)
        psiv_table[i] = ((psiv_table[i - 1] * psiv) % q)

    return [w_table, wv_table, psi_table, psiv_table]

# Merged NTT with pre-processing (optimized) (iterative)
# This is not NTT, this is pre-processing + NTT
# (see: https://eprint.iacr.org/2016/504.pdf)
# A: input polynomial (standard order)
# Psi: 2n-th root of unity
# q: modulus
# B: output polynomial (bit-reversed order)
def NTT(A,Psi,q):
    N = len(A)
    B = [_ for _ in A]

    l = int(math.log(N,2))

    t = N
    m = 1
    while(m<N):
        t = int(t/2)
        for i in range(m):
            j1 = 2*i*t
            j2 = j1 + t - 1
            Psi_pow = intReverse(m+i,l)
            S = pow(Psi,Psi_pow,q)
            for j in range(j1,j2+1):
                U = B[j]
                V = (B[j+t]*S) % q

                B[j]   = (U+V) % q
                B[j+t] = (U-V) % q
        m = 2*m

    return B

# Merged INTT with post-processing (optimized) (iterative)
# This is not NTT, this is pre-processing + NTT
# (see: https://eprint.iacr.org/2016/504.pdf)
# A: input polynomial (Bit-reversed order)
# Psi: 2n-th root of unity
# q: modulus
# B: output polynomial (standard order)
def INTT(A,Psi,q):
    Psi = modinv(Psi,q)
    N = len(A)
    B = [_ for _ in A]

    l = int(math.log(N,2))

    t = 1
    m = N
    while(m>1):
        j1 = 0
        h = int(m/2)
        for i in range(h):
            j2 = j1 + t - 1
            Psi_pow = intReverse(h+i,l)
            S = pow(Psi,Psi_pow,q)
            for j in range(j1,j2+1):
                U = B[j]
                V = B[j+t]

                B[j]   = (U+V) % q
                B[j+t] = (U-V)*S % q
            j1 = j1 + 2*t
        t = 2*t
        m = int(m/2)

    N_inv = modinv(N, q)
    for i in range(N):
        B[i] = (B[i] * N_inv) % q

    return B
