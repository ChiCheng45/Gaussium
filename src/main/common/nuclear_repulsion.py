import numpy as np
from src.main.common import Vector

"""
CoulombsLawMatrix used to calculated the nuclear-nuclear repulsion energies.

DESCRIPTION
    This class contains the methods to calculate the nuclear-nuclear repulsion energies when inputted with a nuclei
    array the method 'create' returns a matrix with all the individual nuclei repulsion energies. The total repulsion
    energy is half the sum of this matrix. E.g.

                HELIUM      HYDROGEN
    HELIUM   [[ 0.          1.36687076]
    HYDROGEN  [ 1.36687076  0.        ]]

ARGUMENTS

    class CoulombsLawMatrix:

        Methods
        -------

        @classmethod
        def create(cls, nuclei_array):

            Parameters
            ----------
            nuclei_array : Nuclei[]
                an array of nuclei objects

            Returns
            -------
            matrix : np.matrix(np.ndarray(int, int), dtype=float)
                a np.matrix contain all individual nuclear-nuclear repulsion energies

        @staticmethod
        def coulombs_law(nuc1, nuc2):

            Parameters
            ----------
            nuc1, nuc2 : Nuclei
                nuclei objects

            Returns
            -------
            ans : float
                the nuclear repulsion energy of nuc1 and nuc2

SEE ALSO
    main.py
    https://en.wikipedia.org/wiki/Coulomb%27s_law

DIAGNOSTICS
    None
"""


class CoulombsLawMatrix:

    @staticmethod
    def coulombs_law(nuc1, nuc2):
        r_12 = Vector.distance(nuc1.coordinates, nuc2.coordinates)
        ans = (nuc1.charge * nuc2.charge) / r_12
        return ans

    @classmethod
    def create(cls, nuclei_array):
        matrix_length = len(nuclei_array)
        matrix = np.matrix(np.zeros((matrix_length, matrix_length)))
        for i in range(matrix_length):
            for j in range(matrix_length):
                if i == j:
                    matrix[i, j] = 0
                else:
                    matrix[i, j] = cls.coulombs_law(nuclei_array[i], nuclei_array[j])
        return matrix


