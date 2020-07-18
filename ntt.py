import math
from helper import *

# --- NTT functions ---

# Iterative NTT (Forward and Inverse)
# arrayIn = input polynomial
# P = modulus
# W = nth root of unity
# inv = 0: forward NTT / 1: inverse NTT
# Input: Standard order/Output: Standard order

def NTT(A, W_table, q):
    # print("DEBUG q:{}".format(q))
    n = len(A)
    B = [x for x in A]

    v = int(math.log(n, 2))

    for i in range(0, v):
        for j in range(0, (2 ** i)):
            for k in range(0, (2 ** (v - i - 1))):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                # w = (W ** ((2 ** i) * k)) % q
                w = W_table[((2 ** i) * k)]

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)

    return B

def INTT(A, W_table, q):
    n = len(A)
    B = [x for x in A]

    v = int(math.log(n, 2))
    
    for i in range(0, v):
        for j in range(0, (2 ** i)):
            for k in range(0, (2 ** (v - i - 1))):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                # w = (W ** ((2 ** i) * k)) % q
                w = W_table[((2 ** i) * k)]

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)

    n_inv = modinv(n, q)
    for i in range(n):
        B[i] = (B[i] * n_inv) % q

    return B

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
