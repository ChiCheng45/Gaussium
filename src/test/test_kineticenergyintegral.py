from unittest import TestCase
from src.main import KineticEnergyIntegral, FileInputNuclei, FileInputBasis
from numpy import testing


class TestKineticEnergyIntegral1(TestCase):

    def setUp(self):
        file_reader_nuclei = FileInputNuclei('HeH+.mol')
        nuclei_array = file_reader_nuclei.create_nuclei_array()
        file_input_basis = FileInputBasis('STO-3G-edited.gbs', nuclei_array)
        basis_set_array = file_input_basis.create_basis_set_array()
        self.kinetic_energy_integral = KineticEnergyIntegral(basis_set_array)

    def test_method_calculate_returns_the_kinetic_energy_of_helium_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 2.1643, 5)

    def test_method_calculate_returs_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 0)
        testing.assert_approx_equal(kinetic_energy, 0.167, 3)

    def test_method_calculate_returs_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 1)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)


class TestKineticEnergyIntegral2(TestCase):

    def setUp(self):
        file_reader_nuclei = FileInputNuclei('HYDROGEN.mol')
        nuclei_array = file_reader_nuclei.create_nuclei_array()
        file_input_basis = FileInputBasis('STO-3G-edited.gbs', nuclei_array)
        basis_set_array = file_input_basis.create_basis_set_array()
        self.kinetic_energy_integral = KineticEnergyIntegral(basis_set_array)

    def test_method_calculate_returns_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)

    def test_method_calculate_returs_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 0)
        testing.assert_approx_equal(kinetic_energy, 0.2365, 4)
