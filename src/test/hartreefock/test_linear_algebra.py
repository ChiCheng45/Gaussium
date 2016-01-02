from src.main.hartreefock import LinearAlgebra
from unittest import TestCase
from numpy import testing
import numpy as np


class TestLinearAlgebra(TestCase):

    def test_transformation_matrix_creates_transformation_matrix(self):
        orbital_overlap = np.matrix([[1.0, 0.4508], [0.4508, 1.0]])
        transformation = LinearAlgebra.transformation_matrix(orbital_overlap)
        testing.assert_array_almost_equal(transformation, np.matrix([[-0.9541, 0.5871], [0.9541, 0.5871]]), decimal=4)

    def test_diagonalize_returns_orbital_energy_eigenvalues(self):
        core_hamiltonian = np.matrix([[], []])
