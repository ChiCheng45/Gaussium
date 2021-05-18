from unittest import TestCase
from unittest.mock import MagicMock, Mock
import numpy as np
from numpy import testing
from gaussium.hartreefock import LinearAlgebra



class TestLinearAlgebra(TestCase):

    def test_transformation_matrix_creates_transformation_matrix(self):
        orbital_overlap = np.array([[1.0, 0.4508], [0.4508, 1.0]])
        linear_algebra = LinearAlgebra(orbital_overlap)

        testing.assert_array_almost_equal(
            linear_algebra.transformation_matrix, np.array([[-0.9541, 0.5871], [0.9541, 0.5871]]), decimal=4
        )

    def test_diagonalize_returns_orbital_energy_eigenvalues_and_orbital_coefficients(self):
        LinearAlgebra.create_transformation_matrix = Mock(return_value=np.array([[-0.9541, 0.5871], [0.9541, 0.5871]]))
        mock_orbital_overlap = MagicMock()
        linear_algebra = LinearAlgebra(mock_orbital_overlap)
        core_hamiltonian = np.array([[-2.6527, -1.3472], [-1.3472, -1.7318]])
        orbital_energies, orbital_coefficients = linear_algebra.diagonalize(core_hamiltonian)

        testing.assert_array_almost_equal(
            orbital_coefficients, np.array([[0.9291, 0.6259], [0.1398, -1.1115]]), decimal=3
        )
        testing.assert_array_almost_equal(orbital_energies, np.array([-2.6742, -1.3043]), decimal=4)
