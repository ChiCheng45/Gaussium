from src.main.objects import RotationSymmetry
from src.main.common import coordinate_distance
import numpy as np
import copy


def check_symmetry(basis1, basis2, basis3, basis4):
    r_1 = basis1.coordinates
    r_2 = basis2.coordinates
    r_3 = basis3.coordinates
    r_4 = basis4.coordinates
    l_1 = basis1.integral_exponents
    l_2 = basis2.integral_exponents
    l_3 = basis3.integral_exponents
    l_4 = basis4.integral_exponents

    if parity(l_1[0]) * parity(l_2[0]) * parity(l_3[0]) * parity(l_4[0]) == -1:
        if r_1[0] == r_2[0] == r_3[0] == r_4[0]:
            return False
        else:
            return True
    elif parity(l_1[1]) * parity(l_2[1]) * parity(l_3[1]) * parity(l_4[1]) == -1:
        if r_1[1] == r_2[1] == r_3[1] == r_4[1]:
            return False
        else:
            return True
    elif parity(l_1[2]) * parity(l_2[2]) * parity(l_3[2]) * parity(l_4[2]) == -1:
        if r_1[2] == r_2[2] == r_3[2] == r_4[2]:
            return False
        else:
            return True
    return True


def parity(num):
    if num % 2 == 0:
        return 1
    else:
        return -1


def sort_index(i, j, k, l):
    a = i
    b = j
    c = k
    d = l
    if a > b:
        a, b = b, a
    if c > d:
        c, d = d, c
    if a > c:
        a, c = c, a
        b, d = d, b
    if a == c:
        if b > d:
            a, c = c, a
            b, d = d, b
    return a, b, c, d


def basis_set_symmetry_matrix(molecule, basis_set):
    symmetry_objects = [None]

    for rotation in molecule.rotation_symmetry:
        rotations = expand_rotation_symmetry(rotation)
        symmetry_objects += rotations
    for reflection in molecule.reflection_symmetry:
        symmetry_objects.append(reflection)

    basis_set_size = len(basis_set)
    symmetry_operation_size = len(symmetry_objects)

    symmetry_matrix = np.empty((basis_set_size + 1, symmetry_operation_size), dtype=object)
    for j in range(len(symmetry_objects)):
        if j == 0:
            symmetry_matrix.itemset((0, j), 'E')
        else:
            symmetry_matrix.itemset((0, j), symmetry_objects[j].symmetry_operation)

    for i in range(1, basis_set_size + 1):
        for j in range(symmetry_operation_size):
            if j == 0:
                symmetry_matrix.itemset((i, j), i)
            else:
                symmetry_matrix.itemset((i, j), symmetry_operation_index(symmetry_objects[j], basis_set,
                basis_set[i - 1]))

    return symmetry_matrix


def expand_rotation_symmetry(rotation_symmetry):
    vector = rotation_symmetry.vector
    fold = rotation_symmetry.fold

    if rotation_symmetry.fold == 2:
        return [rotation_symmetry]

    rotation_operations = []
    for i in range(1, fold):
        rotation_symmetry = RotationSymmetry(fold / i, vector)
        rotation_operations.append(rotation_symmetry)

    return rotation_operations


def symmetry_operation_index(symmetry_operation, basis_set, basis):
    basis_i = copy.deepcopy(basis)
    basis_i.coordinates = symmetry_operation.operate(basis.coordinates)
    x, y, z = symmetry_operation.operate(basis_i.integral_exponents)
    basis_i.integral_exponents = (int(round(x, 1)), int(round(y, 1)), int(round(z, 1)))

    i = 0
    k = 1
    for j, basis_j in enumerate(basis_set):
        if check_basis_functions(basis_i, basis_j):
            k = check_sign(basis_i, basis_j)
            i = j
            break

    return (i + 1) * k


def check_sign(basis_i, basis_x):
    i, j, k = basis_i.integral_exponents
    x, y, z = basis_x.integral_exponents
    i, j, k = even_to_abs(i), even_to_abs(j), even_to_abs(k)

    if i != x or j != y or k != z:
        return -1
    else:
        return 1


def even_to_abs(i):
    if i % 2 == 0:
        return abs(i)
    else:
        return i


def check_basis_functions(basis_i, basis_x):
    i, j, k = basis_i.integral_exponents
    x, y, z = basis_x.integral_exponents
    i, j, k = abs(i), abs(j), abs(k)

    if coordinate_distance(basis_i.coordinates, basis_x.coordinates) <= 1e-3 and i == x and j == y and k == z \
    and check_primitives(basis_i, basis_x):
        return True
    else:
        return False


def check_primitives(basis_i, basis_j):
    primitives_i = basis_i.primitive_gaussian_array
    primitives_j = basis_j.primitive_gaussian_array
    for k in range(len(primitives_i)):
        if primitives_i[k].contraction != primitives_j[k].contraction \
        or primitives_i[k].exponent != primitives_j[k].exponent:
            return False
    return True
