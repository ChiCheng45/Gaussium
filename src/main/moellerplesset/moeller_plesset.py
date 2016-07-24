from src.main.matrixelements import molecular_orbitals
import itertools


class MoellerPlesset:

    def __init__(self, hartree_fock):
        self.hartree_fock = hartree_fock

    def second_order(self):
        electron_energy, orb_energies, orbital_coefficients = self.hartree_fock.begin_scf()

        print('BEGIN MP2 CALCULATION\n')
        occupied_orbitals = self.hartree_fock.electrons // 2
        total_orbitals = len(orb_energies)
        integrals = molecular_orbitals(self.hartree_fock.repulsion, orbital_coefficients)

        correlation = 0
        for i, j in itertools.product(range(occupied_orbitals), repeat=2):
            for a, b in itertools.product(range(occupied_orbitals, total_orbitals), repeat=2):
                out = (orb_energies.item(i) + orb_energies.item(j) - orb_energies.item(a) - orb_energies.item(b))
                correlation += 2 * (integrals.item(i, a, j, b))**2 / out
                correlation -= (integrals.item(i, a, j, b) * integrals.item(i, b, j, a)) / out

        return electron_energy, correlation
