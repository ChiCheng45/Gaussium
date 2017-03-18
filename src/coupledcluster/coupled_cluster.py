import time
from src.common import Indices
from src.coupledcluster import PeturbativeTriples
from src.coupledcluster import SinglesDoubles
from src.matrixelements import molecular_orbitals
from src.matrixelements import spin_basis_anti_physicist
from src.matrixelements import spin_orbital_energies


class CoupledCluster(Indices):

    def __init__(self, hartree_fock, threshold):
        self.hartree_fock_energy, orbital_energies, orbital_coefficients = hartree_fock.begin_scf()
        self.repulsion = spin_basis_anti_physicist(molecular_orbitals(hartree_fock.repulsion, orbital_coefficients))
        self.orbital_energies = spin_orbital_energies(orbital_energies)
        self.occupied_orbitals = hartree_fock.electrons
        self.unoccupied_orbitals = len(self.orbital_energies) - hartree_fock.electrons
        super().__init__(self.occupied_orbitals, self.unoccupied_orbitals)
        self.threshold = threshold


class CoupledClusterSinglesDoubles(CoupledCluster):

    def __init__(self, hartree_fock, threshold=1e-12):
        super().__init__(hartree_fock, threshold)
        self.amplitudes_factory = SinglesDoubles(self.repulsion, self.orbital_energies, self.occupied_orbitals,
        self.unoccupied_orbitals)

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

        return correlation, amplitudes

    def singles_doubles_correlation(self, t):
        correlation = 0
        for i, j, a, b in self.doubles():
            correlation += self.repulsion.item(i, j, a, b) * (0.25 * t[i, j, a, b] + 0.5 * t[i, a] * t[j, b])
        return correlation

    def energies(self):
        correlation = self.calculate_singles_doubles()[0]
        return self.hartree_fock_energy, correlation


class CoupledClusterPerturbativeTriples(CoupledClusterSinglesDoubles):

    def __init__(self, hartree_fock, threshold=1e-12):
        super().__init__(hartree_fock, threshold)
        self.singles_doubles_correlation, self.t_singles_doubles = self.calculate_singles_doubles()
        self.amplitudes_factory = PeturbativeTriples(self.repulsion, self.orbital_energies, self.occupied_orbitals,
        self.unoccupied_orbitals)

    def calculate_perturbative_triples(self):
        print('\n*************************************************************************************************')
        print('\nBEGIN CCSD(T) CALCULATION')
        start = time.clock()
        t_connected, t_disconnected = self.amplitudes_factory.calculate_triples_amplitudes(self.t_singles_doubles)
        correlation_triples = self.perturbative_triples_correlation(t_connected, t_disconnected)
        correlation = self.singles_doubles_correlation + correlation_triples
        print('CCSD(T) SINGLES AND DOUBLES CORRELATION ENERGY: ' + str(self.singles_doubles_correlation) + ' a.u.')
        print('CCSD(T) TRIPLES CORRELATION ENERGY: ' + str(correlation_triples) + ' a.u.')
        print('CCSD(T) CORRELATION ENERGY: ' + str(correlation) + ' a.u.')
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n\n')
        return correlation

    def perturbative_triples_correlation(self, t_connected, t_disconnected):
        correlation_triples = 0
        d = self.amplitudes_factory.denominator
        for i, j, k, a, b, c in self.triples():
            correlation_triples += (1 / 36) * t_connected[i, j, k, a, b, c] * d[i, j, k, a, b, c] \
            * (t_connected[i, j, k, a, b, c] + t_disconnected[i, j, k, a, b, c])
        return correlation_triples

    def energies(self):
        correlation = self.calculate_perturbative_triples()
        return self.hartree_fock_energy, correlation
