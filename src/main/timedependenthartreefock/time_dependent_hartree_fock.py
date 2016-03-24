from src.main.hartreefock import RestrictedHF
from src.main.matrixelements import TDHFMatrix
from src.main.matrixelements import MolecularIntegrals
import numpy as np


class TimeDependentHartreeFock:

    # Under Construction
    def __init__(self):
        pass

    @staticmethod
    def tamm_dancoff_approximation(nuclei_array, basis_set_array, electrons):
        electron_energy, orbital_energies, orbital_coefficients, repulsion = RestrictedHF(nuclei_array, basis_set_array, electrons).begin()

        molecular_integral_matrix = MolecularIntegrals.calculate(repulsion, orbital_coefficients)
        a_matrix = TDHFMatrix(molecular_integral_matrix, orbital_energies, electrons).a_matrix()

        print(a_matrix)

        eigenvalues, unitary = np.linalg.eigh(a_matrix)

        print(eigenvalues)

        return electron_energy
