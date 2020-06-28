import math
from helper import *

# --- NTT functions ---

# Iterative NTT (Forward and Inverse)
# arrayIn = input polynomial
# P = modulus
# W = nth root of unity
# inv = 0: forward NTT / 1: inverse NTT
# Input: Standard order/Output: Standard order

def NTT(A, W, q):
    n = len(A)
    B = [x for x in A]

    v = int(math.log(n, 2))

    for i in range(0, v):
        for j in range(0, (2 ** i)):
            for k in range(0, (2 ** (v - i - 1))):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                w = (W ** ((2 ** i) * k)) % q

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)

    return B

def INTT(A, W, q):
    n = len(A)
    B = [x for x in A]

    v = int(math.log(n, 2))

    for i in range(0, v):
        for j in range(0, (2 ** i)):
            for k in range(0, (2 ** (v - i - 1))):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                w = (W ** ((2 ** i) * k)) % q

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)

    n_inv = modinv(n, q)
    for i in range(n):
        B[i] = (B[i] * n_inv) % q

    return B

"""

# Cooley-Tukey NTT
# arrayIn = input polynomial
# P = modulus
# W = nth root of unity
# inv = 0: forward NTT / 1: inverse NTT
# Input: Standard order/Output: Bit-Reversed order
def CooleyTukeyNTT(arrayIn, P, W, inv):
    N = len(arrayIn)

    if (N == 2):
        arrayOut = [0] * N

        arrayOut[0] = (arrayIn[0] + arrayIn[1]) % P
        arrayOut[1] = (arrayIn[0] - arrayIn[1]) % P

        return arrayOut
    else:
        arrayOut = [0] * N
        w = 1

        arrayIn_even = [0] * (N / 2)
        arrayIn_odd = [0] * (N / 2)

        for i in range(N / 2):
            arrayIn_even[i] = arrayIn[2 * i]
            arrayIn_odd[i] = arrayIn[2 * i + 1]

        arrayOut_even = CooleyTukeyNTT(arrayIn_even, P, (W * W % P), 0)
        arrayOut_odd  = CooleyTukeyNTT(arrayIn_odd, P, (W * W % P), 0)

        for i in range(N / 2):
            arrayOut[i] = (arrayOut_even[i] + w * arrayOut_odd[i]) % P
            arrayOut[i + (N / 2)] = (arrayOut_even[i] - w * arrayOut_odd[i]) % P

            w = w * W

    if (inv):
        N_inv = util.modinv(N, P)
        for i in range(N):
            arrayOut[i] = (arrayOut[i] * N_inv) % P

    return arrayOut

matrix = lambda polynomial, col_length: zip(*[polynomial[i:i + col_length] for i in range(0, len(polynomial), col_length)])

# Four-Step NTT
# arrayIn = input polynomial
# P = modulus
# W = nth root of unity
# inv = 0: forward NTT / 1: inverse NTT
# size = input polynomial partition
# Input: Standard order/Output: Standard order
def FourStepNTT(arrayIn, P, W, inv, size):
    N = len(arrayIn)

    # If this is an inverse transform operation
    if inv:
        N_inv = util.modinv(N, P)
        # Re-order input
        poly = [arrayIn[0]] + list(reversed(arrayIn[1:]))
    else:
        poly = arrayIn

    size0 = size[0]
    size1 = size[1]

    temp0 = 1
    # STEP.1
    if isinstance(size0, list):
        temp0 = util.MultAll(size0)
        STEP_1 = matrix(poly, N/temp0)
        W_0 = (W ** (N/temp0)) % P
        for i in range(N/temp0):
            STEP_1[i] = FourStepNTT(STEP_1[i], P, W_0, 0, size0)
    else:
        temp0 = size0
        STEP_1 = matrix(poly, N/temp0)
        W_0 = (W ** (N/temp0)) % P
        for i in range(N/temp0):
            STEP_1[i] =  CooleyTukeyNTT(STEP_1[i], P, W_0, 0)

    # STEP.2 - Transpose
    STEP_2 = zip(*STEP_1)

    # STEP.3 - Multiply with twiddle factor of N-pt NTT
    STEP_3 = [[0]*(N/temp0)]*size0
    for i in range(temp0):
        STEP_3[i] = [(STEP_2[i][k] * (W ** (i*k)) % P) for k in range(N/temp0)]

    temp1 = 1
    #STEP.4
    if isinstance(size1, list):
        temp1 = util.MultAll(size1)
        W_1 = (W ** (N/temp1)) % P
        for i in range(N/temp1):
            STEP_3[i] = FourStepNTT(STEP_3[i], P, W_1, 0, size1)
    else:
        temp1 = size1
        W_1 = (W ** (N/temp1)) % P
        for i in range(N/temp1):
            STEP_3[i] = CooleyTukeyNTT(STEP_3[i], P, W_1, 0)

    # Final transpose
    STEP_4 = zip(*STEP_3)

    # Convert matrix into array
    STEP_4 = [item for sublist in STEP_4 for item in sublist]

    if inv:
        for i in range(N):
            STEP_4[i] = (STEP_4[i] * N_inv) % P

    return STEP_4
"""
