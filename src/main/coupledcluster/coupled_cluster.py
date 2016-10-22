from src.main.matrixelements import molecular_orbitals
from src.main.matrixelements import spin_basis_anti_physicist
from src.main.matrixelements import spin_orbital_energies
from src.main.coupledcluster import Indices
from src.main.coupledcluster import SinglesDoubles
from src.main.coupledcluster import PeturbativeTriples
import time


class CoupledCluster(Indices):

    def __init__(self, hartree_fock, threshold):
        self.hartree_fock_energy, orbital_energies, orbital_coefficients = hartree_fock.begin_scf()
        self.repulsion = spin_basis_anti_physicist(molecular_orbitals(hartree_fock.repulsion, orbital_coefficients))
        self.orbital_energies = spin_orbital_energies(orbital_energies)
        super().__init__(range(hartree_fock.electrons), range(hartree_fock.electrons, hartree_fock.electrons
        + len(self.orbital_energies) - hartree_fock.electrons))
        self.threshold = threshold


class CoupledClusterSinglesDoubles(CoupledCluster):

    def __init__(self, hartree_fock, threshold=1e-12):
        super().__init__(hartree_fock, threshold)
        self.amplitudes_factory = SinglesDoubles(self.repulsion, self.orbital_energies, self.occupied, self.unoccupied)

    def calculate_singles_doubles(self):
        print('*************************************************************************************************')
        print('\nMP2 INITIAL GUESS')
        start = time.clock()
        amplitudes = self.amplitudes_factory.mp2_initial_guess()
        correlation = self.singles_doubles_correlation(amplitudes)
        print('MP2 CORRELATION ENERGY: ' + str(correlation) + ' a.u.')
        print('TIME TAKEN: ' + str(time.clock() - start) + 's')

        print('\n*************************************************************************************************')
        print('\nBEGIN CCSD ITERATIONS')
        start = time.clock()
        previous_correlation = 0

        while True:

            amplitudes = self.amplitudes_factory.calculate_amplitudes(amplitudes)
            correlation = self.singles_doubles_correlation(amplitudes)
            delta_energy = previous_correlation - correlation
            previous_correlation = correlation
            print('CCSD CORRELATION ENERGY: ' + str(correlation) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n\n')

        return self.hartree_fock_energy, correlation, amplitudes

    def singles_doubles_correlation(self, t):
        correlation = 0
        for i, j, a, b in self.doubles():
            correlation += self.repulsion.item(i, j, a, b) * (0.25 * t[i, j, a, b] + 0.5 * t[i, a] * t[j, b])
        return correlation


class CoupledClusterPerturbativeTriples(CoupledClusterSinglesDoubles):

    def __init__(self, hartree_fock, threshold=1e-12):
        super().__init__(hartree_fock, threshold)
        self.hartree_fock_energy, self.singles_doubles_correlation, self.amplitudes = self.calculate_singles_doubles()
        self.amplitudes_factory = PeturbativeTriples(self.repulsion, self.orbital_energies, self.occupied,
        self.unoccupied)

    def calculate_perturbative_triples(self):
        pass
