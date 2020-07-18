from helper import *
from poly import *


class RNS:
    def __init__(self, n, t, q, b, qnp, bnp):
        self.n = n
        self.q = q  # q = [q0,q1,....]
        self.t = t
        self.b = b  # b = [b0,b1,....]
        self.qnp = qnp
        self.bnp = bnp

        # Constants
        self.m_sk = 2305843009213554689
        self.m_tilde = 1 << 32
        self.gamma = 2305843009213489153

        # Bases
        self.base_q = RNSBase(self.q, qnp)
        self.base_t_gamma = RNSBase([self.t, self.gamma])
        self.base_B = RNSBase(self.b)
        self.base_Bsk = RNSBase(self.b + [self.m_sk])
        self.base_Bsk_m_tilde = RNSBase(self.b + [self.m_sk, self.m_tilde])

        # Base Converters
        self.base_q_to_base_t_gamma = BaseConverter(self.base_q, self.base_t_gamma)
        self.base_q_to_Bsk = BaseConverter(self.base_q, self.base_Bsk)
        self.base_q_to_m_tilde = BaseConverter(self.base_q, RNSBase([self.m_tilde]))
        self.base_B_to_q = BaseConverter(self.base_B, self.base_q)
        self.base_B_to_m_sk = BaseConverter(self.base_B, RNSBase([self.m_sk]))

        ##### Calculate Variables #####

        #### Multiplication Part ####

        # Compute prod(B) mod q
        self.prod_B_mod_q = [self.base_B.base_prod % self.q[i] for i in range(self.base_q.base_size)]

        # Compute prod(q)^(-1) mod Bsk
        self.inv_prod_q_mod_Bsk = [modinv(self.base_q.base_prod % self.base_Bsk[i], self.base_Bsk[i]) for i in
                                   range(self.base_Bsk.base_size)]

        # Compute prod(B) ^ (-1) mod m_sk
        self.inv_prod_B_mod_m_sk = modinv(self.base_B.base_prod % self.m_sk, self.m_sk)

        # Compute m_tilde^(-1) mod Bsk
        self.inv_m_tilde_mod_Bsk = [modinv(self.m_tilde, self.base_Bsk[i]) for i in range(self.base_Bsk.base_size)]

        # Compute prod(q)^(-1) mod m_tilde
        self.inv_prod_q_mod_m_tilde = modinv(self.base_q.base_prod % self.m_tilde, self.m_tilde)

        # Compute prod(q) mod Bsk
        self.prod_q_mod_Bsk = [self.base_q.base_prod % self.base_Bsk[i] for i in range(self.base_Bsk.base_size)]

        #### Decryption part ####

        # Compute gamma^(-1) mod t
        self.inv_gamma_mod_t = modinv(self.gamma % self.t, self.t)

        # Compute prod({t, gamma}) mod q
        self.prod_t_gamma_mod_q = [(self.base_t_gamma[0] * self.base_t_gamma[1] % self.q[i]) for i in
                                   range(self.base_q.base_size)]

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
        # temp = [[0] * self.base_q.base_size for _ in range(self.n)]
        # for i in range(self.n):
        #     for j in range(len(self.q)):
        #         temp[i][j] = (ctext.F[i][j] * self.prod_t_gamma_mod_q[j]) % self.q[j]

        temp = ctext * self.prod_t_gamma_mod_q

        # Convert from q to {t, gamma}
        temp_t_gamma = self.base_q_to_base_t_gamma.fast_convert_array(temp)

        # Multiply by - prod(q) ^ (-1) mod {t, gamma}
        # for i in range(self.n):
        #     for j in range(self.base_t_gamma.base_size):
        #         temp_t_gamma[i][j] = (temp_t_gamma[i][j] * self.neg_inv_q_mod_t_gamma[j]) % self.base_t_gamma[j]
        temp_t_gamma = temp_t_gamma * self.neg_inv_q_mod_t_gamma

        # ptext = [0] * self.n
        ptext = Poly(self.n, self.t)

        # Need to correct values in temp_t_gamma (gamma component only) which are
        # larger than floor(gamma/2)
        gamma_div_2 = self.gamma >> 1
        for i in range(self.n):
            # Need correction because of centered mod
            if temp_t_gamma[1][i] > gamma_div_2:
                ptext[i] = (temp_t_gamma[0][i] + (self.gamma - temp_t_gamma[1][i])) % self.t
            else:
                ptext[i] = (temp_t_gamma[0][i] - temp_t_gamma[1][i]) % self.t
            # If this coefficient was non-zero, multiply by t^(-1)
            # Perform final multiplication by gamma inverse mod t

        ptext = ptext * self.inv_gamma_mod_t

        return ptext

    def divide_and_round_q_last_inplace(self, rns_poly):

        # Add (qi-1)/2 to change from flooring to rounding
        last_modulus = self.base_q[-1]
        half_last_modulus = last_modulus >> 1
        rns_poly[-1] = rns_poly[-1] + half_last_modulus

        for i in range(self.base_q.base_size - 1):
            temp_poly = rns_poly[-1] % self.base_q[i]

            # Change poly modulus for further operations
            temp_poly.q = self.base_q[i]
            half_mod = half_last_modulus % self.base_q[i]
            temp_poly = temp_poly - half_mod

            rns_poly[i] = (rns_poly[i] - temp_poly) * self.inv_q_last_mod_q[i]

        return rns_poly

    """
    # Require: Input in q
    # Ensure: Output in Bsk U {m_tilde}
    """
    def fastbconv_m_tilde(self, rns_poly):

        temp = rns_poly * self.m_tilde

        # new rns poly = [rns_poly in bsk, rns_poly in m_tilde]
        dest = self.base_q_to_Bsk.fast_convert_array(temp)
        dest2 = self.base_q_to_m_tilde.fast_convert_array(temp)

        # Append results in 1 rns_poly
        dest.append(dest2)
        return dest

    """
    # Require: Input should in Bsk U {m_tilde}
    # Ensure: Destination array in Bsk = m U {msk}
    """
    def sm_mrq(self, rns_poly):
        # Compute r_m_tilde
        poly_m_tilde = rns_poly.get_poly(-1)
        poly_m_tilde = -(poly_m_tilde * self.inv_prod_q_mod_m_tilde)

        m_tilde_div_2 = self.m_tilde >> 1;

        for l in range(self.base_Bsk.base_size):
            for i in range(self.n):
                temp = poly_m_tilde[i]
                if temp >= m_tilde_div_2:
                    temp += (self.base_Bsk[l] - self.m_tilde)

                rns_poly[l][i] = ((rns_poly[l][i] + (self.prod_q_mod_Bsk[l] * temp)) * self.inv_m_tilde_mod_Bsk[l]) % \
                                 self.base_Bsk[l]

        rns_poly.drop_last_poly()
        return rns_poly

    """
    # Require: Input1 should in q and Input2 should in Bsk
    # Ensure: Destination array in Bsk = m U {msk}
    """
    def fast_floor(self, rns_poly_q, rns_poly_bsk):

        dest_rns_poly_bsk = self.base_q_to_Bsk.fast_convert_array(rns_poly_q)

        for l in range(self.base_Bsk.base_size):
            for i in range(self.n):
                dest_rns_poly_bsk[l][i] = ((rns_poly_bsk[l][i] + self.base_Bsk[l] - dest_rns_poly_bsk[l][i])
                                           * self.inv_prod_q_mod_Bsk[l]) % self.base_Bsk[l]

        return dest_rns_poly_bsk

    """
    # Require: Input should in Bsk
    # Ensure: Destination array in q
    """
    def fastbconv_sk(self, rns_poly_bsk):
        dest_rns_poly_q = self.base_B_to_q.fast_convert_array(rns_poly_bsk)
        temp_rns_poly_msk = self.base_B_to_m_sk.fast_convert_array(rns_poly_bsk.get_poly(-1))

        alpha_poly_sk = ((temp_rns_poly_msk[0] + self.m_sk) - rns_poly_bsk.get_poly(-1)) * self.inv_prod_B_mod_m_sk

        m_sk_div_2 = self.m_sk >> 1;
        for k in range(self.base_q.base_size):
            for i in range(self.n):
                if alpha_poly_sk[i] > m_sk_div_2:
                    dest_rns_poly_q[k][i] = (dest_rns_poly_q[k][i] +
                                             (self.prod_B_mod_q[k] * (self.m_sk - alpha_poly_sk[i]))) % self.base_q[k]
                else:
                    dest_rns_poly_q[k][i] = (dest_rns_poly_q[k][i] +
                                             (self.m_sk - self.prod_B_mod_q[k] * alpha_poly_sk[i])) % self.base_q[k]

        return dest_rns_poly_q


class BaseConverter:
    def __init__(self, ibase, obase):
        self.ibase = ibase
        self.obase = obase

        self.base_change_matrix = [[0] * obase.base_size for _ in range(ibase.base_size)]
        for i in range(obase.base_size):
            for j in range(ibase.base_size):
                self.base_change_matrix[j][i] = ibase.punctured_prod_array[j] % self.obase[i]

    def fast_convert_array(self, rns_poly):

        # for i in range(self.ibase.base_size):
        #     for k in range(len(rns_poly)):
        #         rns_poly[k][i] = (rns_poly[k][i] * self.ibase.inv_punctured_prod_mod_base_array[i]) % self.ibase[i]
        #
        temp = rns_poly * self.ibase.inv_punctured_prod_mod_base_array

        result_poly = RNSPoly(temp.n, self.obase.rns_base, self.obase.ntt_tables)

        for j in range(self.obase.base_size):
            for k in range(rns_poly.n):
                tmp = 0
                for i in range(self.ibase.base_size):
                    tmp += (temp[i][k] * self.base_change_matrix[i][j])
                result_poly[j][k] = tmp % self.obase[j]
        return result_poly


class RNSBase:
    def __init__(self, base_mod, ntt_tables=None):
        self.rns_base = base_mod.copy()
        self.base_size = len(base_mod)
        if ntt_tables is None:
            self.ntt_tables = [[0, 0, 0, 0] for _ in range(self.base_size)]
        else:
            self.ntt_tables = ntt_tables
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
            temp = ((self.inv_punctured_prod_mod_base_array[i] * rns_value[i]) % self.rns_base[i]) * \
                   self.punctured_prod_array[i]
            value = (value + temp) % self.base_prod
        return value

    """Define a total base component with its sub-modules"""

    def decompose(self, value):
        rns_value = [0] * self.base_size
        for i in range(self.base_size):
            rns_value[i] = value % self.rns_base[i]
        return rns_value
