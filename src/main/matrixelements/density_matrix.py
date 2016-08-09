from numpy import matlib


def density_matrix_restricted(orbital_coefficient, electrons):
    matrix_size = orbital_coefficient.shape[0]
    density_matrix = matlib.zeros((matrix_size, matrix_size))
    for i in range(electrons // 2):
        density_matrix += 2 * orbital_coefficient[:, i] * orbital_coefficient[:, i].T
    return density_matrix


def density_matrix_unrestricted(orbital_coefficient, electrons):
    matrix_size = orbital_coefficient.shape[0]
    density_matrix = matlib.zeros((matrix_size, matrix_size))
    for i in range(electrons):
        density_matrix += orbital_coefficient[:, i] * orbital_coefficient[:, i].T
    return density_matrix


def blocked_density_matrix(orbital_coefficient, electrons_alph, electrons_beta):
    matrix_size = orbital_coefficient.shape[0]
    half_matrix_size = matrix_size // 2

    def density_alph():
        density_matrix_alph = matlib.zeros((matrix_size, matrix_size))
        for i in range(electrons_alph):
            density_matrix_alph += orbital_coefficient[:, i] * orbital_coefficient[:, i].T
        return density_matrix_alph

    def density_beta():
        density_matrix_beta = matlib.zeros((matrix_size, matrix_size))
        for i in range(half_matrix_size, half_matrix_size + electrons_beta):
            density_matrix_beta += orbital_coefficient[:, i] * orbital_coefficient[:, i].T
        return density_matrix_beta

    return density_alph() + density_beta()
