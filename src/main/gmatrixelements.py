from src.main.twoelectronrepulsionintegrals import TwoElectronRepulsion

class TwoElectronPartOfTheFockMatrixElements:

    def __init__(self, density_matrix, basis_set_array):
        self.density_matrix = density_matrix
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        g_ij = 0
        for a in range(0, self.density_matrix.shape[0]):
            for b in range(0, self.density_matrix.shape[0]):
                two_electron_repulsion = TwoElectronRepulsion(self.basis_set_array)
                coulomb_integral = two_electron_repulsion.calculate(i, j, a, b)
                exchange_integral = two_electron_repulsion.calculate(i, b, a, j)
                g_ij += self.density_matrix.item((b, a)) * (coulomb_integral - (1/2) * exchange_integral)
        return g_ij
