from src.main.common import coulomb_matrix
from src.main.common import Symmetry
from src.main.hartreefock import RestrictedHF
from src.main.hartreefock import UnrestrictedHF
from src.main.hartreefock import BlockedHartreeFock
from src.main.kohnsham import RestrictedKohnSham
from src.main.moellerplesset import MoellerPlesset
from src.main.tdhartreefock import TimeDependentHartreeFock
from src.main.coupledcluster import CoupledClusterSinglesDoubles
from src.main.coupledcluster import CoupledClusterPerturbativeTriples
from src.main.objects import PointGroup


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
            ).calculate(tda=False)

        if self.method == 'CIS':
            electron_energy, correlation = TimeDependentHartreeFock(
                RestrictedHF(nuclei_array, basis_set, self.electrons, self.symmetry_object, self.processors)
            ).calculate(tda=True)

        total_energy = electron_energy + nuclear_repulsion + correlation
        print('NUCLEAR REPULSION ENERGY:    ' + str(nuclear_repulsion) + ' a.u.')
        print('SCF ENERGY:                  ' + str(electron_energy) + ' a.u.')
        print('CORRELATION ENERGY:          ' + str(correlation) + ' a.u.')
        print('TOTAL ENERGY:                ' + str(total_energy) + ' a.u.')
        return total_energy
