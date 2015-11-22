from src.main.integrals import NuclearAttractionIntegral


class NuclearAttractionElement:

    def __init__(self, nuclei_array, basis_set_array):
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        v_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
        for a in range(len(primitive_gaussian_array_i)):
            for b in range(len(primitive_gaussian_array_j)):
                c_1 = primitive_gaussian_array_i[a].contraction
                c_2 = primitive_gaussian_array_j[b].contraction
                n_1 = primitive_gaussian_array_i[a].normalisation
                n_2 = primitive_gaussian_array_j[b].normalisation
                for k in range(len(self.nuclei_array)):
                    v_ij += - self.nuclei_array[k].charge * n_1 * n_2 * c_1 * c_2 * NuclearAttractionIntegral.primitive_nuclear_attraction(primitive_gaussian_array_i[a], primitive_gaussian_array_j[b], self.nuclei_array[k])
        return v_ij

