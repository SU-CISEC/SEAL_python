from ntt import *
from generate_prime import *
from random import randint

"""
    # This class never used yet. It may changed in future.
"""
class Modulus:
    def __init__(self, n, mod, root):
        self.mod = mod
        # root = psi
        self.root = root
        self.root_powers = generate_root_powers(n, mod, root)


def generate_modulus(n, t, q_bit):
    # calculate q and psi
    wfound = False
    while (not (wfound)):
        q = generate_large_prime(q_bit)
        # check q = 1 (mod 2n)
        while (not ((q % (2 * n)) == 1)):
            q = generate_large_prime(q_bit)

        # generate NTT parameters
        for i in range(2, q - 1):
            wfound = isrootofunity(i, 2 * n, q)
            if wfound:
                psi = i
                break

    # calculate g
    while (1):
        g = randint(q, 2 ** q_bit - 1)
        if ((g > q) and (gcd(g, q) == 1) and (gcd(g, t) == 1)):
            break

    return q, psi, g