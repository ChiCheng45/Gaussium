from src.main.hartreefock import RestrictedHF
from src.main.moellerplesset import MolecularIntegrals


class MoellerPlesset:

    @staticmethod
    def second_order(nuclei_array, basis_set_array, electrons):
        electron_energy, orbital_energies, orbital_coefficients, repulsion = RestrictedHF(nuclei_array, basis_set_array, electrons).begin()

        correlation = 0
        occupied_orbitals = electrons // 2
        molecular_integral_matrix = MolecularIntegrals.calculate(repulsion, orbital_coefficients)

        print('BEGIN MP2 CALCULATION', end='\n\n')
        for i in range(occupied_orbitals):
            for j in range(occupied_orbitals):
                for a in range(occupied_orbitals, len(orbital_energies)):
                    for b in range(occupied_orbitals, len(orbital_energies)):
                        out1 = (orbital_energies.item(i) + orbital_energies.item(j) - orbital_energies.item(a) - orbital_energies.item(b))
                        out2 = 2 * (molecular_integral_matrix.item(i, a, j, b))**2 / out1
                        out3 = (molecular_integral_matrix.item(i, a, j, b) * molecular_integral_matrix.item(i, b, j, a)) / out1
                        correlation += out2 - out3
        return electron_energy, correlation
