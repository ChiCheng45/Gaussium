from src.main.hartreefock import RestrictedHF
from src.main.moellerplesset import MolecularIntegrals


class MoellerPlesset:

    @staticmethod
    def second_order(nuclei_array, basis_set_array, electrons):
        electron_energy, orbital_energies, orbital_coefficients, repulsion = RestrictedHF(nuclei_array, basis_set_array, electrons).begin()
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
