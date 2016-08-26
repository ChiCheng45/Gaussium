from src.main.matrixelements import molecular_orbitals
from src.main.matrixelements import spin_basis_anti_physicist
from src.main.matrixelements import spin_orbital_energies
from src.main.coupledcluster import SinglesDoublesAmplitudes
import time


class CoupledCluster:

    def __init__(self, hartree_fock, threshold=1e-12):
        self.hartree_fock = hartree_fock
        self.threshold = threshold

    def singles_doubles(self):
        electron_energy, orbital_energies, orbital_coefficients = self.hartree_fock.begin_scf()

        repulsion = spin_basis_anti_physicist(molecular_orbitals(self.hartree_fock.repulsion, orbital_coefficients))
        orbital_energies = spin_orbital_energies(orbital_energies)
        occupied = self.hartree_fock.electrons
        unoccupied = len(orbital_energies) - occupied

        print('*************************************************************************************************')
        print('\nMP2 INITIAL GUESS')
        start = time.clock()
        amplitudes_factory = SinglesDoublesAmplitudes(repulsion, orbital_energies, occupied, unoccupied)
        amplitudes = amplitudes_factory.mp2_initial_guess()
        correlation = amplitudes_factory.calculate_correlation(amplitudes)
        print('MP2 CORRELATION ENERGY: ' + str(correlation) + ' a.u.')
        print('TIME TAKEN: ' + str(time.clock() - start) + 's')

        print('\n*************************************************************************************************')
        print('\nBEGIN CCSD ITERATIONS')
        start = time.clock()
        previous_correlation = 0

        while True:

            amplitudes = amplitudes_factory.calculate_amplitudes(amplitudes)
            correlation = amplitudes_factory.calculate_correlation(amplitudes)
            delta_energy = previous_correlation - correlation
            previous_correlation = correlation
            print('CCSD CORRELATION ENERGY: ' + str(correlation) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n\n')

        return electron_energy, correlation
