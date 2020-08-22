from BFV import *

# Determine n and bit-size of q, then find a
# q satisfying the condition: q = 1 (mod 2n)
#
# Based on n and q, generate NTT parameters
n   = 0
q   = 0
psi = 0

# Determine t
t = 1024

print("--- Starting BFV (SEAL) Demo")

n = 4096

#RNS base q
q = [68719403009, 68719230977, 137438822401]
k = len(q)
q_psi = [24250113, 29008497, 8625844]

# RNS base b = 2 tane rastgele mod sonuncusu m_sk
b = [2305843009213317121, 2305843009213243393, 2305843009213145089]
l = len(b)
b_psi = [829315415491244, 32973993658837, 307554654119321]

# Determine B (bound of distribution X)
sigma = 3.2

# Determine T, p (relinearization)
T = 256
p = 16

# RNS BFV Evaluator
Evaluator = BFV(n, t, q, b, q_psi, b_psi)

# Generate Keys
Evaluator.SecretKeyGen()
Evaluator.PublicKeyGen()

# print secret and public keys
print("--- Secret and Public Key are generated.")
print("* sk: {}".format(Evaluator.sk))
print("* pk[0]: {}".format(Evaluator.pk[0]))
print("* pk[1]: {}".format(Evaluator.pk[1]))
print("")

# print system parameters
print(Evaluator)

# Generate random message
m = 100
# Encode random messages into plaintext polynomials
m_plain, m_len = Evaluator.EncodeInt(m)

# Encrypt message
ct1 = Evaluator.EncryptionSEAL_RNS(m_plain, m_len)
mt = Evaluator.DecryptionSEAL_RNS(ct1)

# Decode decrypted polynomial into integer
nr = Evaluator.DecodeInt(mt, m_len)


if (m_plain == mt) and (nr == m):
    print("--- Encryption-Decryption works.")
else:
    print("--- Encryption-Decryption does not work.")


# Homomorphic Operations
# Evaluator.HomomorphicMultiplication_RNS(ct1, ct1)
#
