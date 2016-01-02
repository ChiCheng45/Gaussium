from src.main.common import Symmetry
from src.main.hartreefock import HartreeFock


class MoellerPlesset:

    @staticmethod
    def second_order(nuclei_array, electrons, multiplicity, basis_set_array):
        electron_energy, orbital_energies, orbital_coefficients, repulsion = HartreeFock(nuclei_array, electrons, multiplicity, basis_set_array).restricted()
        correlation = 0
        occupied_orbitals = electrons // 2
        molecular_integral = MolecularIntegrals(repulsion, orbital_coefficients).calc
        for i in range(occupied_orbitals):
            for j in range(occupied_orbitals):
                for a in range(occupied_orbitals, len(orbital_energies)):
                    for b in range(occupied_orbitals, len(orbital_energies)):
                        out1 = (orbital_energies[i] + orbital_energies[j] - orbital_energies[a] - orbital_energies[b])
                        out2 = 2 * (molecular_integral(i, a, j, b))**2 / out1
                        out3 = (molecular_integral(i, a, j, b) * molecular_integral(i, b, j, a)) / out1
                        correlation += out2 - out3
        return electron_energy, correlation


class MolecularIntegrals:

    def __init__(self, repulsion, orbital_coefficients):
        self.repulsion = repulsion
        self.orbital_coefficients = orbital_coefficients

    def calc(self, i, j, k, l):
        ans = 0
        c = self.orbital_coefficients
        for r in range(c.shape[0]):
            for s in range(c.shape[0]):
                for t in range(c.shape[0]):
                    for u in range(c.shape[0]):
                        coefficients = c.item(r, i) * c.item(s, j) * c.item(t, k) * c.item(u, l)
                        ans += coefficients * self.repulsion[Symmetry.sort(r, s, t, u)]
        return ans
