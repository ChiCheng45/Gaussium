import numpy as np

"""
Matrix

SYNOPSIS
    def __init__(self, matrix_size)
    int matrix_size

    def create_matrix(self, element)
    Element element
    np.matrix matrix

DESCRIPTION
    A class with a create_matrix method that takes in a matrix element calculator and builds the matrix for it for the
    desired matrix_size. For example taking in the orbital overlap matrix element calculator it will build the orbital
    overlap matrix.

ARGUMENTS
    def __init__(self, matrix_size)
    matrix_size Input:  initialise the object with the length of the basis set.
    
    def create_matrix(self, element)
    element     Input:  an object with a calculate method that calculates the matrix element for the position i, j.
    matrix      Output: an numpy matrix, e.g. the overlap matrix.

SEE ALSO
    density_matrix.py
    fock_matrix.py
    kinetic_energy_matrix.py
    nuclear_attraction_matrix.py
    orbital_overlap_matrix.py
    two_electron_repulsion_matrix.py

DIAGNOSTICS
    None, as long as this object and the element objects are correctly initialised.
"""


class Matrix:

    def __init__(self):
        self.matrix_size = 0

    def create_matrix(self, calculate):
        matrix = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
        for i in range(self.matrix_size):
            for j in range(self.matrix_size):
                if i <= j:
                    matrix.itemset((i, j), calculate(i, j))
        return matrix + np.transpose(matrix) - np.diag(np.diag(matrix))
