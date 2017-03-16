from src.common import Indices
from src.matrixelements import molecular_orbitals


class MoellerPlesset(Indices):

    def __init__(self, hartree_fock):
        self.hartree_fock_energy, self.orbital_energies, orbital_coefficients = hartree_fock.begin_scf()
        occupied_orbitals = hartree_fock.electrons // 2
        unoccupied_orbitals = len(self.orbital_energies) - hartree_fock.electrons // 2
        super().__init__(occupied_orbitals, unoccupied_orbitals)
        self.integrals = molecular_orbitals(hartree_fock.repulsion, orbital_coefficients)

    def second_order(self):
        print('BEGIN MP2 CALCULATION\n')
        correlation = 0
        for i, j, a, b in self.restricted_doubles():
            out = self.orbital_energies.item(i) + self.orbital_energies.item(j) - self.orbital_energies.item(a) \
            - self.orbital_energies.item(b)
            correlation += 2 * (self.integrals.item(i, a, j, b))**2 / out
            correlation -= (self.integrals.item(i, a, j, b) * self.integrals.item(i, b, j, a)) / out
        return correlation

    def energies(self):
        correlation = self.second_order()
        return self.hartree_fock_energy, correlation
