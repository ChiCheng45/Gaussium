from src.main.matrixelements import molecular_orbitals
from src.main.matrixelements import spin_basis_anti_physicist
from src.main.matrixelements import spin_orbital_energies
from src.main.coupledcluster import SinglesDoublesAmplitudes
import time


class CoupledCluster:

    def __init__(self, hartree_fock):
        self.hartree_fock = hartree_fock

    def singles_doubles(self):
        electron_energy, orbital_energies, orbital_coefficients = self.hartree_fock.begin_scf()

        repulsion = spin_basis_anti_physicist(molecular_orbitals(self.hartree_fock.repulsion, orbital_coefficients))
        orbital_energies = spin_orbital_energies(orbital_energies)
        occupied_orbitals = self.hartree_fock.electrons
        unoccupied_orbitals = len(orbital_energies) - occupied_orbitals

        amplitudes_factory = SinglesDoublesAmplitudes(repulsion, orbital_energies, occupied_orbitals,
        unoccupied_orbitals)

        amplitudes = amplitudes_factory.mp2_initial_guess()
        correlation = amplitudes_factory.calculate_correlation(amplitudes)

        start = time.clock()
        print('CCSD CORRELATION ENERGY: ' + str(correlation) + ' a.u.')
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n\n')

        return electron_energy, correlation
