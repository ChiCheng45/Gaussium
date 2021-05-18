import os
from unittest import TestCase
from numpy import testing
from gaussium.common import read_mol_file


class TestFileInputNuclei(TestCase):

    def test_create_nuclei_array_returns_array_of_nuclei_with_names(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        nuclei_array, electron_count, multiplicity = read_mol_file(mol_file)
        self.assertEqual(nuclei_array[0].element, 'HELIUM')
        self.assertEqual(nuclei_array[1].element, 'HYDROGEN')

    def test_create_nuclei_array_returns_array_of_nuclei_with_charges(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        nuclei_array, electron_count, multiplicity = read_mol_file(mol_file)
        self.assertEqual(nuclei_array[0].charge, 2)
        self.assertEqual(nuclei_array[1].charge, 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_mass(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        nuclei_array, electron_count, multiplicity = read_mol_file(mol_file)
        self.assertEqual(nuclei_array[0].mass, 4)
        self.assertEqual(nuclei_array[1].mass, 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_coordinates(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        nuclei_array, electron_count, multiplicity = read_mol_file(mol_file)
        testing.assert_array_equal(nuclei_array[0].coordinates, (0.0, 0.0, 0.487732042))
        testing.assert_array_equal(nuclei_array[1].coordinates, (0.0, 0.0, -0.975464083))

    def test_electron_count_return_two(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        nuclei_array, electron_count, multiplicity = read_mol_file(mol_file)
        self.assertEqual(electron_count, 2)
