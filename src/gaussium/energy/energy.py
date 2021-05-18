from gaussium.common import Symmetry
from gaussium.common import coulomb_matrix
from gaussium.coupledcluster import CoupledClusterPerturbativeTriples
from gaussium.coupledcluster import CoupledClusterSinglesDoubles
from gaussium.hartreefock import BlockedHartreeFock
from gaussium.hartreefock import RestrictedHF
from gaussium.hartreefock import UnrestrictedHF
from gaussium.kohnsham import RestrictedKohnSham
from gaussium.moellerplesset import MoellerPlesset
from gaussium.objects import PointGroup
from gaussium.tdhartreefock import TimeDependentHartreeFock
from gaussium.tdhartreefock import TammDancoffApproximation


class Energy:

    def __init__(self, electrons, multiplicity, processors, method):
        self.electrons = electrons
        self.multiplicity = multiplicity
        self.processors = processors
        self.method = method
        self.symmetry_object = Symmetry(PointGroup([], [], [], [], 'C_{1}'), [])

    def calculate_energy(self, nuclei_array, basis_set):
        electron_energy = correlation = 0.0

        coulomb_law_matrix = coulomb_matrix(nuclei_array)
        nuclear_repulsion = coulomb_law_matrix.sum() / 2

        if self.method == 'RHF':
            electron_energy, correlation = RestrictedHF(
                nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors
            ).energies()

        if self.method == 'UHF':
            electron_energy, correlation = UnrestrictedHF(
                nuclei_array, basis_set, self.electrons, self.multiplicity, self.symmetry_object, self.processors
            ).energies()

        if self.method == 'GUHF':
            electron_energy, correlation = BlockedHartreeFock(
                nuclei_array, basis_set, self.electrons, self.multiplicity, self.symmetry_object, self.processors
            ).energies()

        if self.method == 'MP2':
            electron_energy, correlation = MoellerPlesset(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).energies()

        if self.method[0] == 'DFT':
            electron_energy, correlation = RestrictedKohnSham(
                nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors, self.method[1],
                self.method[2]
            ).energies()

        if self.method == 'CCSD':
            electron_energy, correlation = CoupledClusterSinglesDoubles(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).energies()

        if self.method == 'CCSD(T)':
            electron_energy, correlation = CoupledClusterPerturbativeTriples(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).energies()

        if self.method == 'TDHF':
            electron_energy, correlation = TimeDependentHartreeFock(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).calculate()

        if self.method == 'CIS':
            electron_energy, correlation = TammDancoffApproximation(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).calculate()

        total_energy = electron_energy + nuclear_repulsion + correlation
        print('NUCLEAR REPULSION ENERGY:    ' + str(nuclear_repulsion) + ' a.u.')
        print('SCF ENERGY:                  ' + str(electron_energy) + ' a.u.')
        print('CORRELATION ENERGY:          ' + str(correlation) + ' a.u.')
        print('TOTAL ENERGY:                ' + str(total_energy) + ' a.u.')
        return total_energy
