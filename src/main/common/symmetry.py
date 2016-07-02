from src.main.objects import RotationSymmetry
from src.main.common import coordinate_distance
import numpy as np
import copy


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


def non_zero_integral(symmetry_matrix, basis_set, index):
    i, j, k, l = index
    basis_i = basis_set[i]
    basis_j = basis_set[j]
    basis_k = basis_set[k]
    basis_l = basis_set[l]

    symmetry_exponents_i = symmetry_exponents(basis_i.integral_exponents)
    symmetry_exponents_j = symmetry_exponents(basis_j.integral_exponents)
    symmetry_exponents_k = symmetry_exponents(basis_k.integral_exponents)
    symmetry_exponents_l = symmetry_exponents(basis_l.integral_exponents)

    exponents = [symmetry_exponents_i, symmetry_exponents_j, symmetry_exponents_k, symmetry_exponents_l]
    if basis_i.coordinates == basis_j.coordinates == basis_k.coordinates == basis_l.coordinates:

        x = y = z = 0
        for exponent in exponents:
            x += exponent[0]
            y += exponent[1]
            z += exponent[2]

        if x % 2 != 0:
            return False
        if y % 2 != 0:
            return False
        if z % 2 != 0:
            return False
        else:
            return True

    if symmetry_matrix is not None:

        i += 1
        j += 1
        k += 1
        l += 1

        symmetry_operations = symmetry_matrix.shape[1]
        for m in range(1, symmetry_operations):
            a = symmetry_matrix.item(i, m)
            b = symmetry_matrix.item(j, m)
            c = symmetry_matrix.item(k, m)
            d = symmetry_matrix.item(l, m)

            if a > 0 and b > 0 and ((a == i and b == j) or (a == j and b == i)):
                if c < 0 and d < 0 and not phase(symmetry_exponents_k, symmetry_exponents_l) and ((c == -k and d == -l)
                or (c == -l and d == -k)):
                    return False
                elif c < 0 < d == l or k == c > 0 > d:
                    return False
            if c > 0 and d > 0 and ((c == k and d == l) or (c == l and d == k)):
                if a < 0 and b < 0 and not phase(symmetry_exponents_i, symmetry_exponents_j) and ((a == -i and b == -j)
                or (a == -j and b == -i)):
                    return False
                elif a < 0 < b == j or i == a > 0 > b:
                    return False
            if a < 0 and b < 0 and phase(symmetry_exponents_i, symmetry_exponents_j):
                if c < 0 < d == l or k == c > 0 > d:
                    return False
                elif c < 0 and d < 0 and not phase(symmetry_exponents_l, symmetry_exponents_k):
                    return False
            if c < 0 and d < 0 and phase(symmetry_exponents_k, symmetry_exponents_l):
                if a < 0 < b == j or i == a > 0 > b:
                    return False
                elif a < 0 and b < 0 and not phase(symmetry_exponents_i, symmetry_exponents_j):
                    return False
            if a < 0 < b == j and c < 0 < d == l:
                if not phase(symmetry_exponents_i, symmetry_exponents_k):
                    return False
            if i == a > 0 > b and c < 0 < d == l:
                if not phase(symmetry_exponents_j, symmetry_exponents_k):
                    return False
            if a < 0 < b == j and k == c > 0 > d:
                if not phase(symmetry_exponents_i, symmetry_exponents_l):
                    return False
            if i == a > 0 > b and k == c > 0 > d:
                if not phase(symmetry_exponents_j, symmetry_exponents_l):
                    return False

    return True


def phase(exponents_i, exponents_j):

    if exponents_i[0] == exponents_j[0] != 0:
        return True
    if exponents_i[1] == exponents_j[1] != 0:
        return True
    if exponents_i[2] == exponents_j[2] != 0:
        return True

    return False


def symmetry_exponents(integral_exponents):
    i = j = k = 1
    if integral_exponents[0] % 2 == 0:
        i = 0
    if integral_exponents[1] % 2 == 0:
        j = 0
    if integral_exponents[2] % 2 == 0:
        k = 0
    return i, j, k


def basis_set_symmetry_matrix(molecule, basis_set):

    if molecule.point_group == 'C_{2v}' or 'D_{2h}':
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

    else:
        return None


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
        if j == len(basis_set) - 1:
            return 0

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
