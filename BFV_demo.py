from generate_prime import *
from BFV import *

import random

# Parameter generation (or enter manually)

# Determine n and bit-size of q, then find a
# q satisfying the condition: q = 1 (mod 2n)
#
# Based on n and q, generate NTT parameters
#
# Also, determine b, on the order of k*q^2

"""
Parameters for n=1024, q_bit=20, b_bit=50
n  : 1024
q  : 817153
b  : 553590119577601
t  : 256
B  : 10
qnp: [141376, 291861, 376, 241234]
bnp: [129146617481542, 266370398002610, 350750842, 519197795291993]
"""
if __name__ == "__main__":

    n = 2048
    q_bit = 54
    b_bit = 61

    q = [68719403009, 68719230977, 137438822401]
    qnp = []
    roots = [24250113, 29008497, 8625844]
    for i, psi in enumerate(roots):
        psiv = modinv(psi, q[i])
        w = pow(psi, 2, q[i])
        wv = modinv(w, q[i])
        params = [psi, psiv, w, wv]
        qnp.append(params)

    # b and bnp variables are not correct yet. UPDATE this for MULT.
    b = [2305843009213317121, 2305843009213243393, 2305843009213145089]
    bnp = [129146617481542, 266370398002610, 350750842, 519197795291993]

    # print((q % (2 * n)) == 1)

    # # calculate q and qnp
    # while(1):
    #     # q = generate_large_prime(q_bit)
    #     # # check q = 1 (mod 2n)
    #     # while (not ((q % (2*n)) == 1)):
    #     #     q = generate_large_prime(q_bit)
    #
    #     # generate NTT parameters
    #     for i in range(2,q-1):
    #         if pow(i,2*n,q) == 1:
    #             if pow(i,n,q) == (q-1):
    #                 pru = [i**x % q for x in range(1,2*n)]
    #                 if not(1 in pru):
    #                     qnp[0] = (i**2 % q)
    #                     qnp[1] = modinv(qnp[0],q)
    #                     qnp[2] = i
    #                     qnp[3] = modinv(qnp[2],q)
    #                     break
    #             else:
    #                 continue
    #             break
    #         else:
    #             continue
    #         break
    #     else:
    #         continue
    #     break

    # # calculate b and bnp
    # while(1):
    #     # b = q*generate_large_prime(b_bit-q_bit)
    #     # # check q = 1 (mod 2n)
    #     # while (not ((b % (2*n)) == 1)):
    #     #     b = q*generate_large_prime(b_bit-q_bit)
    #
    #     # generate NTT parameters
    #     for i in range(2,b-1):
    #         if pow(i,2*n,b) == 1:
    #             if pow(i,n,b) == (b-1):
    #                 pru = [i**x % b for x in range(1,2*n)]
    #                 if not(1 in pru):
    #                     bnp[0] = (i**2 % b)
    #                     bnp[1] = modinv(bnp[0],b)
    #                     bnp[2] = i
    #                     bnp[3] = modinv(bnp[2],b)
    #                     break
    #             else:
    #                 continue
    #             break
    #         else:
    #             continue
    #         break
    #     else:
    #         continue
    #     break
    print('qnp = ', qnp)
    print('bnp = ', bnp)

    # Determine t
    t = 1024

    # Determine B (bound of distribution X)
    B = 6

    # Determine T, p (relinearization)
    T = 1024
    p = 16

    # Generate BFV evaluator
    Evaluator = BFV(n, q, t, b, B, qnp, bnp)

    print(Evaluator)

    # Generate Keys
    Evaluator.SecretKeyGen()
    Evaluator.PublicKeyGen()
    # Evaluator.EvalKeyGenV1(T)
    # Evaluator.EvalKeyGenV2(p)

    # Generate random message
    # m = Poly(n,t,np)
    # m.randomize(2)
    # mr1 = random.randint(0,t-1)
    # mr2 = random.randint(0,t-1)
    m = 10
    # mr2 = 20

    # Encode message
    # m1  = Evaluator.Encode(mr1)
    # m2  = Evaluator.Encode(mr2)

    # Encrypt message
    ct = Evaluator.Encryption(m)

    # Add two message
    # ct = Evaluator.HomomorphicAddition(ct1, ct2)

    # Multiply two message
    # ct = Evaluator.HomomorphicMultiplication(ct1,ct2)
    # ct = Evaluator.RelinearizationV1(ct)

    # Decrypt message
    r = Evaluator.Decryption(ct)
    # r = Evaluator.DecryptionV2(ct)

    print(r)

    # Decode message
    rm = Evaluator.Decode(r)

    # if(mr1*mr2) == rm:
    #     print("Homomorphic Multiplication is working correctly.")

    #TEST MULT
    print(rm)
    print(mr1)
    print(mr2)


    # MULT is not working. Check it!!!




    #
