import numpy as np


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
        : np.matrix

        """
        print(type(function))
        matrix = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
        for i in range(self.matrix_size):
            for j in range(self.matrix_size):
                if i <= j:
                    matrix.itemset((i, j), function(i, j))
        return matrix + np.transpose(matrix) - np.diag(np.diag(matrix))
