from generate_prime import *
from BFV import *
from helper import *

from random import randint

# Parameter generation (or enter manually)

# Determine n and bit-size of q, then find a
# q satisfying the condition: q = 1 (mod 2n)
#
# Based on n and q, generate NTT parameters

# Parameter (generated or pre-defined)
PC = 1 # 0: generate -- 1: pre-defined

n   = 0
q   = 0
psi = 0
psiv= 0
w   = 0
wv  = 0

g   = 0
# g (gamma) parameter is a constant 61-bit integer in SEAL
# For a K-bit q, gamma could be a K-bit integer larger than q and co-prime to q, t

# Determine t
t = 1024

print("--- Starting BFV (SEAL) Demo")

if PC:
    n = 4096

    q_bit = 37

    # Parameters for q
    # q   = 137438691329
    # psi = 22157790
    # psiv= modinv(psi,q)
    # w   = pow(psi,2,q)
    # wv  = modinv(w,q)

    #RNS base q
    q = [68719403009, 68719230977, 137438822401]
    k = len(q)
    q_psi = [24250113, 29008497, 8625844]
    q_psiv = [modinv(q_psi[i], q[i]) for i in range(k)]
    q_w = [pow(q_psi[i], 2, q[i]) for i in range(k)]
    q_wv = [modinv(q_w[i], q[i]) for i in range(k)]

    # RNS base b = 2 tane rastgele mod sonuncusu m_sk
    b = [2305843009213317121, 2305843009213243393, 2305843009213145089]
    l = len(b)
    b_psi = [829315415491244, 32973993658837, 307554654119321]
    b_psiv = [modinv(b_psi[i], b[i]) for i in range(l)]
    b_w = [pow(b_psi[i], 2, b[i]) for i in range(l)]
    b_wv = [modinv(b_w[i], b[i]) for i in range(l)]

    g   = 0x1fffffffffc80001 # original integer from SEAL

    print("* q is calculated.")
    print("* n   : {}".format(n))
    print("* q   : {}".format(q))
    print("* Parameters are calculated.")
    print("* w   : {}".format(w))
    print("* wv  : {}".format(wv))
    print("* psi : {}".format(psi))
    print("* psiv: {}".format(psiv))
    print("* g   : {}".format(g))

else:
    n = 4096

    q_bit = 25

    # calculate q and qnp
    wfound = False
    while(not(wfound)):
        q = generate_large_prime(q_bit)
        # check q = 1 (mod 2n)
        while (not ((q % (2*n)) == 1)):
            q = generate_large_prime(q_bit)

        # generate NTT parameters
        for i in range(2,q-1):
            wfound = isrootofunity(i,2*n,q)
            if wfound:
                psi = i
                psiv= modinv(psi,q)
                w   = pow(psi,2,q)
                wv  = modinv(w,q)
                break

    # calculate g
    while(1):
        g = randint(q,2**q_bit - 1)
        if((g>q) and (gcd(g,q) == 1) and (gcd(g,t) == 1)):
            break

    print("* q is calculated.")
    print("* n   : {}".format(n))
    print("* q   : {}".format(q))
    print("* Parameters are calculated.")
    print("* w   : {}".format(w))
    print("* wv  : {}".format(wv))
    print("* psi : {}".format(psi))
    print("* psiv: {}".format(psiv))
    print("* g   : {}".format(g))

# Determine B (bound of distribution X)
sigma = 3.1
B = int(10*sigma)

# Determine T, p (relinearization)
T = 256
p = 16

# # Generate tables
# w_table    = [1]*n
# wv_table   = [1]*n
# psi_table  = [1]*n
# psiv_table = [1]*n
# for i in range(1,n):
#     w_table[i]    = ((w_table[i-1]   *w)    % q)
#     wv_table[i]   = ((wv_table[i-1]  *wv)   % q)
#     psi_table[i]  = ((psi_table[i-1] *psi)  % q)
#     psiv_table[i] = ((psiv_table[i-1]*psiv) % q)
#
# qnp = [w_table,wv_table,psi_table,psiv_table]

#Generate RNS NTT tables
q_ntt_tables = []
for j in range(k):
    w_table    = [1]*n
    wv_table   = [1]*n
    psi_table  = [1]*n
    psiv_table = [1]*n
    for i in range(1,n):
        w_table[i]    = ((w_table[i-1]   *q_w[j])    % q[j])
        wv_table[i]   = ((wv_table[i-1]  *q_wv[j])   % q[j])
        psi_table[i]  = ((psi_table[i-1] *q_psi[j])  % q[j])
        psiv_table[i] = ((psiv_table[i-1]*q_psiv[j]) % q[j])

    qnp = [w_table,wv_table,psi_table,psiv_table]
    q_ntt_tables.append(qnp)

b_ntt_tables = []
for j in range(l):
    w_table    = [1]*n
    wv_table   = [1]*n
    psi_table  = [1]*n
    psiv_table = [1]*n
    for i in range(1,n):
        w_table[i]    = ((w_table[i-1]   *b_w[j])    % b[j])
        wv_table[i]   = ((wv_table[i-1]  *b_wv[j])   % b[j])
        psi_table[i]  = ((psi_table[i-1] *b_psi[j])  % b[j])
        psiv_table[i] = ((psiv_table[i-1]*b_psiv[j]) % b[j])
    bnp = [w_table, wv_table, psi_table, psiv_table]
    b_ntt_tables.append(bnp)

# Generate BFV evaluator
# Evaluator = BFV(n, q, t, B, qnp)
# RNS BFV Evaluator
Evaluator = BFV(n, t, B, q, q_ntt_tables, b, b_ntt_tables)

# Generate Keys
Evaluator.SecretKeyGen()
Evaluator.PublicKeyGen()
# Evaluator.EvalKeyGenV1(T)
# Evaluator.EvalKeyGenV2(p)

# print secret and public keys
print("--- Secret and Public Key are generated.")
print("* sk: {}".format(Evaluator.sk))
print("* pk[0]: {}".format(Evaluator.pk[0]))
print("* pk[1]: {}".format(Evaluator.pk[1]))
print("")

# print system parameters
print(Evaluator)

# Generate random message
n1 = randint(0, Evaluator.q_total-1)

print("--- Random integer n1 is generated.")
print("* n1: {}".format(n1))
print("")

# Encode random messages into plaintext polynomials
m_plain, m_value = Evaluator.EncodeSEAL(31)

print("--- n1 is encoded as polynomial m1(x).")
print("* m1(x): {}".format(m_plain))
print("* m1_value: {}".format(m_value))
print("")

# Encrypt message
ct1 = Evaluator.EncryptionSEAL_RNS(m_plain, m_value)

print("--- m1 is encrypted as ct1.")
print("* ct1[0]: {}".format(ct1[0]))
print("* ct1[1]: {}".format(ct1[1]))
print("")

# Decrypt message
# mt = Evaluator.DecryptionSEAL(ct1,g)
mt = Evaluator.DecryptionSEAL_RNS(ct1)

print("--- ct1 is decrypted as mt(x).")
print("* mt(x): {}".format(mt))
print("")

# Decode decrypted polynomial into integer
nr = Evaluator.DecodeSEAL(mt)

print("--- mt(x) is decoded as integer nr.")
print("* nr: {}".format(nr))
print("")

if (m_plain == mt) and (nr == n1):
    print("--- Encryption-Decryption works.")
else:
    print("--- Encryption-Decryption does not work.")


# Homomorphic Operations
# Evaluator.HomomorphicMultiplication_RNS(ct1, ct1)
#
