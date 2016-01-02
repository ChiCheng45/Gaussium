from src.main.hartreefock import LinearAlgebra
from unittest import TestCase
from numpy import testing
import numpy as np


class TestLinearAlgebra(TestCase):

    def test_transformation_matrix_creates_transformation_matrix(self):
        orbital_overlap = np.matrix([[1.0, 0.4508], [0.4508, 1.0]])
        transformation = LinearAlgebra.transformation_matrix(orbital_overlap)
        testing.assert_array_almost_equal(transformation, np.matrix([[-0.9541, 0.5871], [0.9541, 0.5871]]), decimal=4)

    def test_diagonalize_returns_orbital_energy_eigenvalues_and_orbital_coefficients(self):
        transformation = np.matrix([[-0.9541, 0.5871], [0.9541, 0.5871]])
        core_hamiltonian = np.matrix([[-2.6527, -1.3472], [-1.3472, -1.7318]])
        orbital_energies, orbital_coefficients = LinearAlgebra(transformation).diagonalize(core_hamiltonian)
        testing.assert_array_almost_equal(orbital_coefficients, np.matrix([[0.9291, 0.6259], [0.1398, -1.1115]]), decimal=3)
        testing.assert_array_almost_equal(orbital_energies, np.array([-2.6742, -1.3043]), decimal=4)
