from unittest import TestCase
from src.main import FileInputBasis
from src.main import FileInputNuclei


class TestFileInputBasis(TestCase):

    def setUp(self):
        self.file_reader_nuclei = FileInputNuclei('HeH+.mol')
        self.nuclei_array = self.file_reader_nuclei.create_nuclei_array()

    def test_create_basis_set_array_returns_basis_array_with_correct_names_for_sto3g(self):
        file_reader_basis = FileInputBasis('STO-3G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_name(),'HELIUM')
        self.assertEquals(basis_array[1].get_name(),'HYDROGEN')

    def test_create_basis_set_array_returns_basis_array_with_correct_coordinates_for_sto3g(self):
        file_reader_basis = FileInputBasis('STO-3G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_x(), 0)
        self.assertEquals(basis_array[1].get_x(), 0)
        self.assertEquals(basis_array[0].get_y(), 0)
        self.assertEquals(basis_array[1].get_y(), 0)
        self.assertEquals(basis_array[0].get_z(), 0.7313240)
        self.assertEquals(basis_array[1].get_z(), -0.7313240)

    def test_create_basis_set_array_returns_basis_array_with_orbital_type_for_sto3g(self):
        file_reader_basis = FileInputBasis('STO-3G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_orbital_type(), 'S')
        self.assertEquals(basis_array[1].get_orbital_type(), 'S')

    def test_create_basis_set_array_returns_basis_array_with_array_of_coefficients_for_sto3g(self):
        file_reader_basis = FileInputBasis('STO-3G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_array_of_coefficients(), [[0.15432897, 6.36242139], [0.53532814, 1.15892300], [0.44463454, 0.31364979]])
        self.assertEquals(basis_array[1].get_array_of_coefficients(), [[0.15432897, 3.42525091], [0.53532814, 0.62391373], [0.44463454, 0.16885540]])

    def test_create_basis_set_array_returns_basis_array_with_correct_names_for_321g(self):
        file_reader_basis = FileInputBasis('3-21G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_name(),'HELIUM')
        self.assertEquals(basis_array[1].get_name(),'HELIUM')
        self.assertEquals(basis_array[2].get_name(),'HYDROGEN')
        self.assertEquals(basis_array[3].get_name(),'HYDROGEN')


    def test_create_basis_set_array_returns_basis_array_with_correct_coordinates_for_321g(self):
        file_reader_basis = FileInputBasis('3-21G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()

        self.assertEquals(basis_array[0].get_x(), 0)
        self.assertEquals(basis_array[1].get_x(), 0)
        self.assertEquals(basis_array[2].get_x(), 0)
        self.assertEquals(basis_array[3].get_x(), 0)

        self.assertEquals(basis_array[0].get_y(), 0)
        self.assertEquals(basis_array[1].get_y(), 0)
        self.assertEquals(basis_array[2].get_y(), 0)
        self.assertEquals(basis_array[3].get_y(), 0)

        self.assertEquals(basis_array[0].get_z(), 0.7313240)
        self.assertEquals(basis_array[1].get_z(), 0.7313240)
        self.assertEquals(basis_array[2].get_z(), -0.7313240)
        self.assertEquals(basis_array[3].get_z(), -0.7313240)

    def test_create_basis_set_array_returns_basis_array_with_orbital_type_for_321g(self):
        file_reader_basis = FileInputBasis('3-21G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_orbital_type(), 'S')
        self.assertEquals(basis_array[1].get_orbital_type(), 'S')
        self.assertEquals(basis_array[2].get_orbital_type(), 'S')
        self.assertEquals(basis_array[3].get_orbital_type(), 'S')

    def test_create_basis_set_array_returns_basis_array_with_array_of_coefficients_for_sto3g(self):
        file_reader_basis = FileInputBasis('3-21G.gbs', self.nuclei_array)
        basis_array = file_reader_basis.create_basis_set_array()
        self.assertEquals(basis_array[0].get_array_of_coefficients(), [[0.1752300, 13.6267000], [0.8934830, 1.9993500]])
        self.assertEquals(basis_array[1].get_array_of_coefficients(), [[1.0000000, 0.3829930]])
        self.assertEquals(basis_array[2].get_array_of_coefficients(), [[0.1562850, 5.4471780], [0.9046910, 0.8245470]])
        self.assertEquals(basis_array[3].get_array_of_coefficients(), [[1.0000000, 0.1831920]])