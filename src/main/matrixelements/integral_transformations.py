import numpy as np


def molecular_orbitals(repulsion, coefficients):
    matrix_size = coefficients.shape[0]
    for r in range(matrix_size):
        for s in range(matrix_size):
            repulsion[r, s, :, :] = np.transpose(coefficients) * repulsion[r, s, :, :] * coefficients
    for t in range(matrix_size):
        for u in range(matrix_size):
            repulsion[:, :, t, u] = np.transpose(coefficients) * repulsion[:, :, t, u] * coefficients
    return repulsion


def spin_basis_set(repulsion):
    half_matrix_size = repulsion.shape[0]
    matrix_size = half_matrix_size * 2
    zero_matrix = np.zeros((half_matrix_size, half_matrix_size))
    spin_orbital_repulsion = np.zeros((matrix_size, matrix_size, matrix_size, matrix_size))
    for r in range(matrix_size):
        for s in range(matrix_size):

            if r < half_matrix_size and s < half_matrix_size:
                block_matrix = np.bmat([
                        [repulsion[r, s, :, :], zero_matrix],
                        [zero_matrix, repulsion[r, s, :, :]]])
                spin_orbital_repulsion[r, s, :, :] = block_matrix
            elif r >= half_matrix_size and s >= half_matrix_size:
                block_matrix = np.bmat([
                        [repulsion[r - half_matrix_size, s - half_matrix_size, :, :], zero_matrix],
                        [zero_matrix, repulsion[r - half_matrix_size, s - half_matrix_size, :, :]]])
                spin_orbital_repulsion[r, s, :, :] = block_matrix
            else:
                block_matrix = np.bmat([
                        [zero_matrix, zero_matrix],
                        [zero_matrix, zero_matrix]])
                spin_orbital_repulsion[r, s, :, :] = block_matrix

    return spin_orbital_repulsion
