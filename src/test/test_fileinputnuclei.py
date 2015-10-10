from unittest import TestCase
from src.main import FileInputNuclei


class TestFileInputBasis(TestCase):

    def setUp(self):
        self.file_reader_nuclei = FileInputNuclei('HeH+.mol')


    def test_create_nuclei_array_returns_array_of_nuclei_with_names(self):
        nuclei_array = self.file_reader_nuclei.create_nuclei_array()
        self.assertEquals(nuclei_array[0].get_name(), 'HELIUM')
        self.assertEquals(nuclei_array[1].get_name(), 'HYDROGEN')

    def test_create_nuclei_array_returns_array_of_nuclei_with_charges(self):
        nuclei_array = self.file_reader_nuclei.create_nuclei_array()
        self.assertEquals(nuclei_array[0].get_charge(), 2)
        self.assertEquals(nuclei_array[1].get_charge(), 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_mass(self):
        nuclei_array = self.file_reader_nuclei.create_nuclei_array()
        self.assertEquals(nuclei_array[0].get_mass(), 4)
        self.assertEquals(nuclei_array[1].get_mass(), 1)

    def test_create_nuclei_array_returns_array_of_nuclei_with_coordinates(self):
        nuclei_array = self.file_reader_nuclei.create_nuclei_array()
        self.assertEquals(nuclei_array[0].get_x(), 0)
        self.assertEquals(nuclei_array[1].get_x(), 0)
        self.assertEquals(nuclei_array[0].get_y(), 0)
        self.assertEquals(nuclei_array[1].get_y(), 0)
        self.assertEquals(nuclei_array[0].get_z(), 0.7313240)
        self.assertEquals(nuclei_array[1].get_z(), -0.7313240)