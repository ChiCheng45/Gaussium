from src.main.matrixelements import molecular_orbitals


class MoellerPlesset:

    def __init__(self, hartree_fock):
        self.hartree_fock = hartree_fock
        self.electrons = self.hartree_fock.electrons

    def second_order(self):
        electron_energy, orbital_energies, orbital_coefficients, repulsion = self.hartree_fock.begin_scf()

        correlation = 0
        occupied_orbitals = self.electrons // 2
        molecular_integral = molecular_orbitals(repulsion, orbital_coefficients)

        print('BEGIN MP2 CALCULATION', end='\n\n')

        for i in range(occupied_orbitals):
            for j in range(occupied_orbitals):
                for a in range(occupied_orbitals, len(orbital_energies)):
                    for b in range(occupied_orbitals, len(orbital_energies)):
                        out1 = (orbital_energies.item(i) + orbital_energies.item(j) - orbital_energies.item(a)
                        - orbital_energies.item(b))
                        out2 = 2 * (molecular_integral.item(i, a, j, b))**2 / out1
                        out3 = (molecular_integral.item(i, a, j, b) * molecular_integral.item(i, b, j, a)) / out1
                        correlation += out2 - out3

        return electron_energy, correlation
