import numpy as np


class DIIS:
    """Direct Inversion in the Iterative Subspace for improved SCF convergence.

    Attributes
    ----------
    matrix_size : int
    overlap : np.matrix
    linear_algebra : LinearAlgebra
        Object containing the transformation matrix and linear algebra related methods
    fock_array : List
        A list of guess fock matrices.
    error_array : List
        A list of the error vectors.
    diis_subspace : int
        Max size of the DIIS subspace can get before it gets cleared and restarted.
    begin : bool
        Check to turn DIIS on/off depending on the error vectors.

    References
    ----------
    P. Pulay, J. Comput. Chem. 3, 556 (1982).

    """
    def __init__(self, overlap, linear_algebra, diis_subspace=8):
        self.matrix_size = overlap.shape[0]
        self.overlap = overlap
        self.linear_algebra = linear_algebra
        self.fock_array = []
        self.error_array = []
        self.diis_subspace = diis_subspace
        self.begin = False

    def fock_matrix(self, fock, density):
        """Runs DIIS and return a the Fock matrix or a DIIS optimized Fock matrix.

        Parameters
        ----------
        fock : np.matrix
        density : np.matrix

        Returns
        -------
        fock : np.matrix

        """
        error = self.overlap * density * fock - fock * density * self.overlap
        error = self.linear_algebra.orthonormalize(error)

        if np.any(error < 0.1 * np.ones((self.matrix_size, self.matrix_size))):
            self.begin = True
        if np.all(error < 1e-6 * np.ones((self.matrix_size, self.matrix_size))):
            self.begin = False

        if self.begin:
            self.error_array.append(error)
            self.fock_array.append(fock)

            if len(self.fock_array) > self.diis_subspace:  # DIIS subspace reset
                self.fock_array = []
                self.error_array = []
                return fock

            if len(self.fock_array) > 1:
                fock = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
                diis_coefficients = self.create_b_matrix()

                for l in range(len(self.fock_array)):
                    fock += diis_coefficients[l, 0] * self.fock_array[l]

        return fock

    def create_b_matrix(self):
        """Creates the B Matrix for DIIS optimization.

        Occasionally there can be an issues with the b_matrix becoming singularThis methods corrects itself by
        removing a fock matrix and error vector and trying again.

        Returns
        -------
        diis_coefficients : np.matrix

        """
        array_length = len(self.fock_array)

        b_matrix = np.matrix(np.zeros((array_length + 1, array_length + 1)))
        vector_k = np.matrix(np.zeros(array_length + 1))
        vector_k.itemset(array_length, 1)

        for i in range(array_length):
            b_matrix.itemset((i, array_length), 1)
            b_matrix.itemset((array_length, i), 1)
            for j in range(array_length):
                error_i = self.error_array[i]
                error_j = self.error_array[j]
                b_matrix.itemset((i, j), np.trace(error_i * np.transpose(error_j)))

        try:
            diis_coefficients = np.linalg.solve(b_matrix, np.transpose(vector_k))
            return diis_coefficients
        except np.linalg.linalg.LinAlgError:
            self.fock_array.pop(0)
            self.error_array.pop(0)
            print("np.linalg.linalg.LinAlgError: removing a DIIS vector")
            return self.create_b_matrix()
