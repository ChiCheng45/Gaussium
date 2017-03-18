from unittest import TestCase
from numpy import testing
from src.common import read_mol_file


class TestFileInputNuclei(TestCase):

    def test_create_nuclei_array_returns_array_of_nuclei_with_names(self):
        nuclei_array, electron_count, multiplicity = read_mol_file('HeH+.mol')
        self.assertEquals(nuclei_array[0].element, 'HELIUM')
        self.assertEquals(nuclei_array[1].element, 'HYDROGEN')

    def test_create_nuclei_array_returns_array_of_nuclei_with_charges(self):
        nuclei_array, electron_count, multiplicity = read_mol_file('HeH+.mol')
        self.assertEquals(nuclei_array[0].charge, 2)
        self.assertEquals(nuclei_array[1].charge, 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_mass(self):
        nuclei_array, electron_count, multiplicity = read_mol_file('HeH+.mol')
        self.assertEquals(nuclei_array[0].mass, 4)
        self.assertEquals(nuclei_array[1].mass, 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_coordinates(self):
        nuclei_array, electron_count, multiplicity = read_mol_file('HeH+.mol')
        testing.assert_array_equal(nuclei_array[0].coordinates, (0.0, 0.0, 0.487732042))
        testing.assert_array_equal(nuclei_array[1].coordinates, (0.0, 0.0, -0.975464083))

    def test_electron_count_return_two(self):
        nuclei_array, electron_count, multiplicity = read_mol_file('HeH+.mol')
        self.assertEquals(electron_count, 2)
