import numpy as np
import itertools


class Matrix:
    """Class for creating matrices efficiently by only calculating the triangular matrix and reflecting.

    Attributes
    ----------
    matrix_size : int

    """

    def __init__(self, matrix_size):
        self.matrix_size = matrix_size

    def create_matrix(self, function):
        """Create matrix method that takes in a function to calculate the matrix elements.

        Parameters
        ----------
        function : method

        Returns
        -------
        : np.array

        """
        matrix = np.zeros((self.matrix_size, self.matrix_size))
        for i, j in itertools.product(range(self.matrix_size), repeat=2):
            if i <= j:
                matrix.itemset((i, j), function(i, j))
        return matrix + matrix.T - np.diag(np.diag(matrix))
