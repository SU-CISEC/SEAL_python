from helper import *
from poly import *


class RNS:
    def __init__(self, n, q, t, b, qnp, bnp):
        self.n = n
        self.q = q  # q = [q0,q1,....]
        self.qnp = qnp
        self.t = t
        self.b = b
        self.bnp = bnp

        # Constants
        self.m_sk = 2305843009213554689
        self.m_tilde = 1 << 32
        self.gamma = 2305843009213489153

        # Bases
        self.base_q = RNSBase(self.q)
        self.base_t_gamma = RNSBase([self.t, self.gamma])
        self.base_B = RNSBase(self.b)
        self.base_Bsk = RNSBase(self.b + [self.m_sk])
        self.base_Bsk_m_tilde = RNSBase(self.b + [self.m_sk, self.m_tilde])

        # Base Converters
        self.base_q_to_base_t_gamma = BaseConverter(self.base_q, self.base_t_gamma)
        self.base_q_to_Bsk = BaseConverter(self.base_q, self.base_Bsk)
        self.base_q_to_m_tilde = BaseConverter(self.base_q, RNSBase([self.m_tilde]))
        self.base_B_to_q = BaseConverter(self.base_B, self.base_q)
        self.base_B_to_m_sk = BaseConverter(self.base_B,  RNSBase([self.m_sk]))

        ##### Calculate Variables #####

        #### Multiplication Part ####

        # Compute prod(B) mod q
        self.prod_B_mod_q = [0] * self.base_q.base_size
        for i in range(self.base_q.base_size):
            self.prod_B_mod_q[i] = self.base_B.base_prod % self.q[i]

        # Compute prod(q)^(-1) mod Bsk
        self.inv_prod_q_mod_Bsk = [0] * self.base_Bsk.base_size
        for i in range(self.base_Bsk.base_size):
            self.inv_prod_q_mod_Bsk[i] = self.base_q.base_prod % self.base_Bsk[i]
            self.inv_prod_q_mod_Bsk[i] = modinv(self.inv_prod_q_mod_Bsk[i], self.base_Bsk[i])

        # Compute prod(B) ^ (-1) mod m_sk
        self.inv_prod_B_mod_m_sk = self.base_B.base_prod % self.m_sk
        self.inv_prod_B_mod_m_sk = modinv(self.inv_prod_B_mod_m_sk, self.m_sk)

        # Compute m_tilde^(-1) mod Bsk
        self.inv_m_tilde_mod_Bsk = [0] * self.base_Bsk.base_size
        for i in range(self.base_Bsk.base_size):
            self.inv_m_tilde_mod_Bsk[i] = modinv(self.m_tilde, self.base_Bsk[i])

        # Compute prod(q)^(-1) mod m_tilde
        self.inv_prod_q_mod_m_tilde = self.base_q.base_prod % self.m_tilde
        self.inv_prod_q_mod_m_tilde = modinv(self.inv_prod_q_mod_m_tilde, self.m_tilde)

        #### Decryption part ####

        # Compute gamma^(-1) mod t
        self.inv_gamma_mod_t = modinv(self.gamma % self.t, self.t)

        # Compute prod({t, gamma}) mod q
        self.prod_t_gamma_mod_q = [0] * self.base_q.base_size
        for i in range(self.base_q.base_size):
            self.prod_t_gamma_mod_q[i] = self.base_t_gamma[0] * self.base_t_gamma[1] % self.q[i]

        # Compute -prod(q)^(-1) mod {t, gamma}
        self.neg_inv_q_mod_t_gamma = [0] * self.base_t_gamma.base_size
        for i in range(self.base_t_gamma.base_size):
            self.neg_inv_q_mod_t_gamma[i] = self.base_q.base_prod % self.base_t_gamma[i]
            self.neg_inv_q_mod_t_gamma[i] = self.base_t_gamma[i] - modinv(self.neg_inv_q_mod_t_gamma[i],
                                                                          self.base_t_gamma[i])

        # Compute q[last]^(-1) mod q[i] for i = 0..last-1
        # This is used by modulus switching and rescaling
        self.inv_q_last_mod_q = [0] * (self.base_q.base_size - 1)
        for i in range(self.base_q.base_size - 1):
            self.inv_q_last_mod_q[i] = modinv(self.base_q[self.base_q.base_size - 1], self.base_q[i])

    def decrypt_scale_and_round(self, ctext):

        # Compute | gamma * t | _qi * ct(s)
        temp = [[0] * self.base_q.base_size for _ in range(self.n)]
        for i in range(self.n):
            for j in range(len(self.q)):
                temp[i][j] = (ctext.F[i][j] * self.prod_t_gamma_mod_q[j]) % self.q[j]

        # Convert from q to {t, gamma}
        temp_t_gamma = self.base_q_to_base_t_gamma.fast_convert_array(temp)

        # Multiply by - prod(q) ^ (-1) mod {t, gamma}
        for i in range(self.n):
            for j in range(self.base_t_gamma.base_size):
                temp_t_gamma[i][j] = (temp_t_gamma[i][j] * self.neg_inv_q_mod_t_gamma[j]) % self.base_t_gamma[j]

        ptext = [0] * self.n

        # Need to correct values in temp_t_gamma (gamma component only) which are
        # larger than floor(gamma/2)
        gamma_div_2 = self.gamma >> 1
        for i in range(self.n):
            # Need correction because of centered mod
            if temp_t_gamma[i][1] > gamma_div_2:
                ptext[i] = (temp_t_gamma[i][0] + (self.gamma - temp_t_gamma[i][1])) % self.t
            else:
                ptext[i] = (temp_t_gamma[i][0] - temp_t_gamma[i][1]) % self.t
            # If this coefficient was non-zero, multiply by t^(-1)
            # Perform final multiplication by gamma inverse mod t
            ptext[i] = (ptext[i] * self.inv_gamma_mod_t) % self.t

        return ptext
class BaseConverter:
    def __init__(self, ibase, obase):
        self.ibase = ibase
        self.obase = obase

        self.base_change_matrix = [[0] * obase.base_size for _ in range(ibase.base_size)]
        for i in range(obase.base_size):
            for j in range(ibase.base_size):
                self.base_change_matrix[j][i] = ibase.punctured_prod_array[j] % self.obase[i]

    def fast_convert_array(self, poly):

        for i in range(self.ibase.base_size):
            for k in range(len(poly)):
                poly[k][i] = (poly[k][i] * self.ibase.inv_punctured_prod_mod_base_array[i]) % self.ibase[i]

        result_poly = [[0] * self.obase.base_size for _ in range(len(poly))]

        for j in range(self.obase.base_size):
            for k in range(len(poly)):
                tmp = 0
                for i in range(self.ibase.base_size):
                    tmp += (poly[k][i] * self.base_change_matrix[i][j])
                result_poly[k][j] = tmp % self.obase[j]
        return result_poly


class RNSBase:
    def __init__(self, base_mod):
        self.rns_base = base_mod.copy()
        self.base_size = len(base_mod)
        # Initialize base variables
        self.base_prod = 1
        self.punctured_prod_array = [1] * self.base_size
        self.inv_punctured_prod_mod_base_array = [1] * self.base_size

        # Calculate base variables
        if self.base_size > 1:
            # Calculate base product
            for i in range(self.base_size):
                self.base_prod *= self.rns_base[i]

            # Calculate punctured product array
            for i in range(self.base_size):
                temp = 1
                for j in range(self.base_size):
                    if i != j:
                        temp *= self.rns_base[j]
                self.punctured_prod_array[i] = temp

            # Calculate inverse punctured product array
            for i in range(self.base_size):
                self.inv_punctured_prod_mod_base_array[i] = modinv(self.punctured_prod_array[i], self.rns_base[i])

        else:
            self.base_prod = self.rns_base[0]
            self.punctured_prod_array[0] = 1
            self.inv_punctured_prod_mod_base_array[0] = 1

    def __getitem__(self, index):
        return self.rns_base[index]

    def compose(self, rns_value):
        value = 0
        for i in range(self.base_size):
            temp = ((self.inv_punctured_prod_mod_base_array[i] * rns_value[i]) % self.rns_base[i]) * self.punctured_prod_array[i]
            value = (value + temp) % self.base_prod
        return value

    """Define a total base component with its sub-modules"""
    def decompose(self, value):
        rns_value = [0] * self.base_size
        for i in range(self.base_size):
            rns_value[i] = value % self.rns_base[i]
        return rns_value
