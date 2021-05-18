import numpy as np
from gaussium.common import coordinate_distance


def coulomb_matrix(nuclei_array):
    """Creates a matrix containing the repulsion energies of each
    nuclei-nuclei interaction.

    Parameters
    ----------
    nuclei_array : List[Nuclei]

    Returns
    -------
    matrix : np.array

    """
    matrix_length = len(nuclei_array)
    matrix = np.zeros((matrix_length, matrix_length))
    for i in range(matrix_length):
        for j in range(matrix_length):
            if i != j:
                matrix[i, j] = coulombs_law(nuclei_array[i], nuclei_array[j])
    return matrix


def coulombs_law(nuc1, nuc2):
    """Computes the nuclear-nuclear repulsion energy for two nuclei.

    Parameters
    ----------
    nuc1 : Nuclei
    nuc2 : Nuclei

    Returns
    -------
    ans : float

    """
    r_12 = coordinate_distance(nuc1.coordinates, nuc2.coordinates)
    ans = (nuc1.charge * nuc2.charge) / r_12
    return ans
