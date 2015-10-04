from unittest import TestCase
from src.main import FileInputBasis
from src.main import FileInputNuclei


class TestFileInputBasis(TestCase):

    def setUp(self):
        self.file_reader_basis = FileInputBasis('3-21G.gbs')
        self.file_reader_nuclei = FileInputNuclei('HeH+.mol')
        self.nuclei_array = self.file_reader_nuclei.create_nuclei_array()

    def test_create_basis_set_array_returns(self):
        basis_array = self.file_reader_basis.create_basis_set_array_2(self.nuclei_array)
        i = 0
        print(basis_array[i].get_name())
        print(basis_array[i].get_array_of_coefficients())

    # def test_create_basis_set_calls_len_method(self):
    #     basis_array = self.file_reader_basis.create_basis_set_array_2(self.nuclei_array)