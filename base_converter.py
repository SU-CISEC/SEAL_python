# We will first generate all required parameters for base convertion.
# For hardware, we will load all parameters from here
# For software, we will load parameters from memory

class BaseConverter:
    def __init__(self):

        ############Constant Parameters##########

        # n = coeff_count = poly_modulus_degree = 4096
        self.n = 4096
        # t = plain_modulus = 1024
        self.t = 1024

        self.m_sk = 0x1fffffffffe00001
        self.m_tilda = 1 << 32
        self.gamma = 0x1fffffffffc80001

        # Array of coefficient small moduli
        # q = coeff_modulus = [68719403009, 68719230977]
        self.q_arr = [68719403009, 68719230977]

        # Array of auxiliary moduli
        self.aux_base_array_ = [0x1fffffffffb40001, 0x1fffffffff500001]

        # Array of auxiliary U {m_sk_} moduli
        # b = bsk_modulus = new base modulus
        self.bsk_modulus = [0x1fffffffffb40001, 0x1fffffffff500001, 0x1fffffffffe00001]

        # Array of plain modulus U gamma
        self.plain_gamma_array_ = [self.t, self.gamma]

        ############Calculated Parameters##########

        # Punctured products of the coeff moduli
        self.coeff_products_array_ = [68719230977, 68719403009]

        # Matrix which contains the products of coeff moduli mod aux
        self.coeff_base_products_mod_aux_bsk_array_ = [[68719230977, 68719403009], [68719230977, 68719403009], [68719230977, 68719403009]]

        # Array of inverse coeff modulus products mod each small coeff mods
        self.inv_coeff_base_products_mod_coeff_array_ = [26179219651, 42540076863]

        # Array of coeff moduli products mod m_tilde
        self.coeff_base_products_mod_mtilde_array_ = [4294721537, 4294893569]

        # Array of coeff modulus products times m_tilda mod each coeff modulus
        self.mtilde_inv_coeff_base_products_mod_coeff_array_ = [42540557849, 26178779624]

        # Matrix of the inversion of coeff modulus products mod each auxiliary mods
        self.inv_coeff_products_all_mod_aux_bsk_array_ = [1779959502326169890, 245485185625809909]

        # Matrix of auxiliary mods products mod each coeff modulus
        self.aux_base_products_mod_coeff_array_ = [68676968414, 68683522014, 68703633290, 68710186890]

        # Array of inverse auxiliary mod products mod each auxiliary mods
        self.inv_aux_base_products_mod_aux_array_ = [737870114790509117, 1567972894413747653]

        # Array of auxiliary bases products mod m_sk_
        self.aux_base_products_mod_msk_array_ = [2305843009202159617, 2305843009208713217]

        # Coeff moduli products inverse mod m_tilde
        self.inv_coeff_products_mod_mtilde_ = 2349129729

        # Auxiliary base products mod m_sk_ (m1*m2*...*ml)-1 mod m_sk
        self.inv_aux_products_mod_msk_ = 954944906924050867

        # # Gamma inverse mod plain modulus
        # self.inv_gamma_mod_plain_ = 0

        # Auxiliary base products mod coeff moduli (m1*m2*...*ml) mod qi
        self.aux_products_all_mod_coeff_array_ = [48397954621, 54976261965]

        # Array of m_tilde inverse mod Bsk = m U {msk}
        self.inv_mtilde_mod_bsk_array_ = [2303168996393096849, 2299650559177685249, 2304717108767884289]

        # Array of all coeff base products mod Bsk
        self.coeff_products_all_mod_bsk_array_ = [2283888126783854594, 2283888140199073794, 2283888120881158146]

        # # Matrix of coeff base product mod plain modulus and gamma
        # self.coeff_products_mod_plain_gamma_array_
        #
        # # Array of negative inverse all coeff base product mod plain modulus and gamma
        # self.neg_inv_coeff_products_all_mod_plain_gamma_array_
        #
        # # Array of plain_gamma_product mod coeff base moduli
        # self.plain_gamma_product_mod_coeff_array_
        #
        # # For modulus switching: inverses of the last coeff base modulus
        # self.inv_last_coeff_mod_array_

    def fastbconv_mtilde(self, encrypted):
        # [[0]*k]*n
        temp_coeff_transition = [[0] * len(self.q_arr)] * (self.n)
        # Compute in Bsk first; we compute |m_tilde*q^-1i| mod qi
        for i, q_i in enumerate(self.q_arr):
            for c, enc in enumerate(encrypted):
                temp_coeff_transition[i][c] = (enc[i] * self.mtilde_inv_coeff_base_products_mod_coeff_array_[i]) % q_i

        # [[0]*(l+1)]*n
        destination = [[0] * (len(self.bsk_modulus) + 1)] * (self.n)
        for j, b_j in enumerate(self.bsk_modulus):
            for c, enc in enumerate(encrypted):
                res_sum = 0
                for i, q_i in enumerate(self.q_arr):
                    # Product is 60 bit + 61 bit = 121 bit, so can sum up to 127 of them with no reduction
                    # Thus need coeff_base_mod_count_ <= 127
                    # Reduction is needed
                    temp = temp_coeff_transition[i][c] * coeff_base_products_mod_aux_bsk_array_[j][i]
                    res_sum += temp
                destination[j, c] = res_sum % b_j

        # Computing the last element (mod m_tilde) and add it at the end of destination array
        for c in range(self.n):
            res_sum = 0
            for i in range(len(self.q_arr)):
                # Product is 60 bit + 33 bit = 93 bit
                res_sum += (temp_coeff_transition[i][c] * self.coeff_base_products_mod_mtilde_array_[i])

            destination[len(self.bsk_modulus), c] = res_sum % self.m_tilda

        return destination

    """
        Require: Input should in Bsk U {m_tilde}
        Ensure: Destination array in Bsk = m U {msk}
    """

    def mont_rq(self, in_dest):
        # [[0]*l]*n
        l = len(self.bsk_modulus)
        destination = [[0] * l] * (self.n)

        for j, b_j in enumerate(self.bsk_modulus):
            for c in range(n):
                r_mtilde = (in_dest[l][c] * self.inv_coeff_products_mod_mtilde_) % self.m_tilda

                temp = self.coeff_products_all_mod_bsk_array_[j] * r_mtilde
                temp += in_dest[j][c]
                r_mtilde = temp % b_j
                destination[j][c] = (self.inv_mtilde_mod_bsk_array_[j] * r_mtilde) % b_j

    """
        Require: Input in q U m U {msk}
        Ensure: Destination array in Bsk
    """

    def fast_floor(self, in_dest):
        dest = self.fastbconv(in_dest)
        l = len(self.bsk_modulus)
        destination = [[0] * l] * (self.n)

        for i, b_i in enumerate(self.bsk_modulus):
            for c in range(self.n):
                destination[i][c] = (in_dest[i][c] + b_i - dest[i][c]) * self.inv_coeff_products_all_mod_aux_bsk_array_[
                    i]

        return destination

    """
        Require: Input in base Bsk = M U {msk}
        Ensure: Output in base q
    """

    def fastbconv_sk(self, in_dest):
        a = len(self.aux_base_array_)

        # [[0]*a]*n
        temp_coeff_transition = [[0] * a] * (self.n)
        for i, a_i in enumerate(self.aux_base_array_):
            for c in range(self.n):
                temp_coeff_transition[i][c] = (in_dest[i][c] * self.inv_aux_base_products_mod_aux_array_[i]) % a_i

        # [[0]*(k)]*n
        destination = [[0] * (len(self.coeff_products_array_))] * (self.n)
        for j, q_j in enumerate(self.coeff_products_array_):
            for c, enc in enumerate(in_dest):
                res_sum = 0
                for i, a_i in enumerate(self.aux_base_array_):
                    # Product is 60 bit + 61 bit = 121 bit, so can sum up to 127 of them with no reduction
                    # Thus need coeff_base_mod_count_ <= 127
                    # Reduction is needed
                    temp = temp_coeff_transition[i][c] * self.aux_base_products_mod_coeff_array_[j][i]
                    res_sum += temp
                destination[j][c] = res_sum % q_j


        # Compute alpha_sk
        # Require: Input is in Bsk
        # we only use coefficient in B
        # Fast convert B -> m_sk
        tmp = [0] * self.n
        for c in range(self.n):
            res_sum = 0
            for i, a_i in enumerate(self.aux_base_array_):
                # Product is 60 bit + 61 bit = 121 bit, so can sum up to 127 of them with no reduction
                # Thus need coeff_base_mod_count_ <= 127
                # Reduction is needed
                temp = temp_coeff_transition[i][c] * self.aux_base_products_mod_msk_array_[i]
                res_sum += temp
            tmp[c] = res_sum % self.m_sk

        alpha_sk = [0] * self.n
        # x_sk is allocated in input[aux_base_mod_count_]
        for c in range(self.n):
            negated_input = self.m_sk - in_dest[a][c]
            alpha_sk[c] = ((tmp[c] + negated_input) * self.inv_aux_products_mod_msk_) % self.m_sk


        m_sk_div_2 = self.m_sk >> 1
        for i, q_i in enumerate(self.q_arr):
            for c in range(self.n):

                res_sum = 0
                # Correcting alpha_sk since it is a centered modulo
                if (alpha_sk[c] > m_sk_div_2):
                    res_sum = aux_products_all_mod_coeff_array_ [i] * (self.m_sk - alpha_sk[c])
                    destination[i][c] =  destination[i][c] + res_sum % q_i
                # No correction needed
                else:
                    res_sum = (q_i - aux_products_all_mod_coeff_array_[i]) * alpha_sk[c]
                    destination[i][c] = destination[i][c] + res_sum % q_i

        return destination

    '''
        Require: Input in q
        Ensure: Output in Bsk = {m1,...,ml} U {msk}
    '''

    def fastbconv(self, in_dest):
        # [[0]*k]*n
        temp_coeff_transition = [[0] * len(self.q_arr)] * (self.n)
        # Compute in Bsk first; we compute |m_tilde*q^-1i| mod qi
        for i, q_i in enumerate(self.q_arr):
            for c, enc in enumerate(encrypted):
                temp_coeff_transition[i][c] = (enc[i] * self.mtilde_inv_coeff_base_products_mod_coeff_array_[i]) % q_i

        # [[0]*(l)]*n
        destination = [[0] * (len(self.bsk_modulus))] * (self.n)
        for j, b_j in enumerate(self.bsk_modulus):
            for c, enc in enumerate(encrypted):
                res_sum = 0
                for i, q_i in enumerate(self.q_arr):
                    # Product is 60 bit + 61 bit = 121 bit, so can sum up to 127 of them with no reduction
                    # Thus need coeff_base_mod_count_ <= 127
                    # Reduction is needed
                    temp = temp_coeff_transition[i][c] * coeff_base_products_mod_aux_bsk_array_[j][i]
                    res_sum += temp
                destination[j][c] = res_sum % b_j

        return destination

    def print_params(self):
        pass
