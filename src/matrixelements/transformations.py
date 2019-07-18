import numpy as np
import itertools


def molecular_orbitals(repulsion, coefficients):
    """Converts the two electron repulsion integrals from atomic orbital to the molecular orbital basis set.

    Parameters
    ----------
    repulsion : np.array
    coefficients : np.matrix
        This is the orbital coefficients matrix from the Hartree Fock equations.

    Returns
    -------
    repulsion : np.array

    """
    matrix_size = coefficients.shape[0]
    for r, s in itertools.product(range(matrix_size), repeat=2):
        repulsion[r, s, :, :] = np.transpose(coefficients) @ repulsion[r, s, :, :] @ coefficients
    for t, u in itertools.product(range(matrix_size), repeat=2):
        repulsion[:, :, t, u] = np.transpose(coefficients) @ repulsion[:, :, t, u] @ coefficients
    return repulsion


def blocked_spin_basis_set(repulsion):
    """Converts the two electron repusion integrals to spatial orbitals to spin orbitals in the blocked form.

    Parameters
    ----------
    repulsion : np.array

    Returns
    -------
    repulsion : np.array

    """
    half_matrix_size = repulsion.shape[0]
    matrix_size = half_matrix_size * 2
    zero_matrix = np.zeros((half_matrix_size, half_matrix_size))
    spin_orbital_repulsion = np.zeros((matrix_size, matrix_size, matrix_size, matrix_size))

    for r, s in itertools.product(range(matrix_size), repeat=2):

        # r == alph spin and s == alph spin
        if r < half_matrix_size and s < half_matrix_size:
            spin_orbital_repulsion[r, s, :, :] = np.block([
                    [repulsion[r, s, :, :], zero_matrix],
                    [zero_matrix, repulsion[r, s, :, :]]
            ])

        # r == beta spin and s == beta spin
        elif r >= half_matrix_size and s >= half_matrix_size:
            spin_orbital_repulsion[r, s, :, :] = np.block([
                    [repulsion[r - half_matrix_size, s - half_matrix_size, :, :], zero_matrix],
                    [zero_matrix, repulsion[r - half_matrix_size, s - half_matrix_size, :, :]]
            ])

        # r == alph and s == beta or r == beta and s == alph
        else:
            spin_orbital_repulsion[r, s, :, :] = np.block([
                    [zero_matrix, zero_matrix],
                    [zero_matrix, zero_matrix]
            ])

    return spin_orbital_repulsion


def spin_basis_set(repulsion):
    """Converts the two electron repusion integrals of spatial orbitals to spin orbitals.

    Parameters
    ----------
    repulsion : np.array

    Returns
    -------
    repulsion : np.array

    """
    matrix_size = repulsion.shape[0] * 2
    spin_orbital_repulsion = np.zeros((matrix_size, matrix_size, matrix_size, matrix_size))

    for r, s, t, u in itertools.product(range(matrix_size), repeat=4):
        if r % 2 == s % 2 and t % 2 == u % 2:
            spin_orbital_repulsion.itemset((r, s, t, u), repulsion.item(r // 2, s // 2, t // 2, u // 2))

    return spin_orbital_repulsion


def spin_orbital_energies(orbital_energies):
    """Converts the the orbital energies for spatial orbital to spin orbitals.

    Parameters
    ----------
    orbital_energies : List

    Returns
    -------
    orbital_energies : List

    """
    energy_spin_orbital = []
    for i in range(len(orbital_energies)):
        for j in range(2):
            energy_spin_orbital.append(orbital_energies[i])
    return energy_spin_orbital


def spin_basis_anti_physicist(repulsion):
    """Converts the two electron repusion integrals of spatial orbitals to antisymmetrized two electron integrals of
    spin orbitals in the physicists notation.

    Parameters
    ----------
    repulsion : np.array

    Returns
    -------
    spin_orbital_repulsion : np.array

    """
    matrix_size = repulsion.shape[0] * 2
    spin_orbital_repulsion = np.zeros((matrix_size, matrix_size, matrix_size, matrix_size))

    for r, s, t, u in itertools.product(range(matrix_size), repeat=4):
        out1 = (r % 2 == t % 2) * (s % 2 == u % 2) * repulsion.item(r // 2, t // 2, s // 2, u // 2)
        out2 = (r % 2 == u % 2) * (s % 2 == t % 2) * repulsion.item(r // 2, u // 2, s // 2, t // 2)
        spin_orbital_repulsion.itemset((r, s, t, u), out1 - out2)

    return spin_orbital_repulsion
