from src.main.objects import RotationSymmetry
from src.main.common import coordinate_distance
from src.main.common import vector_add
import numpy as np
import copy


class Symmetry:

    def __init__(self, molecule, basis_set_array):
        self.point_group = molecule.point_group
        self.basis_set = basis_set_array
        self.symmetry_objects = self.symmetry_object_list()
        self.int_operate_dict = self.symmetry_exponent_dict()
        self.symmetry_matrix = self.basis_set_symmetry_matrix()

    def none_zero_integral(self, index):
        i, j, k, l = index
        basis_i = self.basis_set[i]
        basis_j = self.basis_set[j]
        basis_k = self.basis_set[k]
        basis_l = self.basis_set[l]

        exponents_kl = vector_add(basis_k.integral_exponents, basis_l.integral_exponents)
        exponents_ij = vector_add(basis_i.integral_exponents, basis_j.integral_exponents)
        symmetry_exponents_ijkl = self.symmetry_exponents(vector_add(exponents_ij, exponents_kl))

        if basis_i.coordinates == basis_j.coordinates == basis_k.coordinates == basis_l.coordinates:
            if symmetry_exponents_ijkl.count(0) != 3:
                return False
            else:
                return True

        i += 1
        j += 1
        k += 1
        l += 1

        symmetry_operations = self.symmetry_matrix.shape[1]
        for m in range(1, symmetry_operations):
            a = self.symmetry_matrix.item(i, m)
            b = self.symmetry_matrix.item(j, m)
            c = self.symmetry_matrix.item(k, m)
            d = self.symmetry_matrix.item(l, m)

            w = abs(a)
            x = abs(b)
            y = abs(c)
            z = abs(d)

            if (a, b, c, d).count(0) == 0 and ((i == w and j == x) or (i == x and j == w)) and ((k == y and l == z)
            or (k == z and l == y)) and symmetry_exponents_ijkl != self.int_operate_dict[(m, symmetry_exponents_ijkl)]:
                return False

        return True

    def symmetry_exponent_dict(self):
        operate_dict = {}
        for i in range(1, len(self.symmetry_objects)):
            operate_dict[(i, (0, 0, 0))] = (0, 0, 0)
            operate_dict[(i, (1, 0, 0))] = self.symmetry_objects[i].int_operate((1, 0, 0))
            operate_dict[(i, (0, 1, 0))] = self.symmetry_objects[i].int_operate((0, 1, 0))
            operate_dict[(i, (0, 0, 1))] = self.symmetry_objects[i].int_operate((0, 0, 1))
            operate_dict[(i, (1, 1, 0))] = self.symmetry_objects[i].int_operate((1, 1, 0))
            operate_dict[(i, (1, 0, 1))] = self.symmetry_objects[i].int_operate((1, 0, 1))
            operate_dict[(i, (0, 1, 1))] = self.symmetry_objects[i].int_operate((0, 1, 1))
            operate_dict[(i, (1, 1, 1))] = self.symmetry_objects[i].int_operate((1, 1, 1))
        return operate_dict

    def symmetry_exponents(self, integral_exponents):
        i = j = k = 1
        if integral_exponents[0] % 2 == 0:
            i = 0
        if integral_exponents[1] % 2 == 0:
            j = 0
        if integral_exponents[2] % 2 == 0:
            k = 0
        return i, j, k

    def symmetry_object_list(self):
        symmetry_objects = [None]
        for rotation in self.point_group.rotation_symmetry:
            rotations = self.expand_rotation_symmetry(rotation)
            symmetry_objects += rotations
        for reflection in self.point_group.reflection_symmetry:
            symmetry_objects.append(reflection)
        symmetry_objects += self.point_group.inversion_symmetry
        return symmetry_objects

    def basis_set_symmetry_matrix(self):

        basis_set_size = len(self.basis_set)
        symmetry_operation_size = len(self.symmetry_objects)

        symmetry_matrix = np.empty((basis_set_size + 1, symmetry_operation_size), dtype=object)
        for j in range(len(self.symmetry_objects)):
            if j == 0:
                symmetry_matrix.itemset((0, j), 'E')
            else:
                symmetry_matrix.itemset((0, j), self.symmetry_objects[j].symmetry_operation)

        for i in range(1, basis_set_size + 1):
            for j in range(symmetry_operation_size):
                if j == 0:
                    symmetry_matrix.itemset((i, j), i)
                else:
                    symmetry_matrix.itemset((i, j), self.symmetry_operation_index(self.symmetry_objects[j],
                    self.basis_set[i - 1]))

        return symmetry_matrix

    def expand_rotation_symmetry(self, rotation_symmetry):
        vector = rotation_symmetry.vector
        fold = rotation_symmetry.fold

        if rotation_symmetry.fold == 2:
            return [rotation_symmetry]

        rotation_operations = []
        for i in range(1, fold):
            rotation_symmetry = RotationSymmetry(fold / i, vector)
            rotation_operations.append(rotation_symmetry)

        return rotation_operations

    def symmetry_operation_index(self, symmetry_operation, basis):
        basis_i = copy.deepcopy(basis)
        basis_i.coordinates = symmetry_operation.operate(basis.coordinates)
        basis_i.integral_exponents = symmetry_operation.int_operate(basis_i.integral_exponents)

        i = k = 0
        for j, basis_j in enumerate(self.basis_set):
            if self.check_basis_functions(basis_i, basis_j):
                k = self.check_sign(basis_i, basis_j)
                i = j
                break
            if j == len(self.basis_set) - 1:
                return 0

        return (i + 1) * k

    def check_sign(self, basis_i, basis_x):
        i, j, k = basis_i.integral_exponents
        x, y, z = basis_x.integral_exponents
        i, j, k = self.even_to_abs(i), self.even_to_abs(j), self.even_to_abs(k)

        if i != x or j != y or k != z:
            return -1
        else:
            return 1

    def even_to_abs(self, i):
        if i % 2 == 0:
            return abs(i)
        else:
            return i

    def check_basis_functions(self, basis_i, basis_x):
        i, j, k = basis_i.integral_exponents
        x, y, z = basis_x.integral_exponents
        i, j, k = abs(i), abs(j), abs(k)

        if coordinate_distance(basis_i.coordinates, basis_x.coordinates) <= 1e-3 and i == x and j == y and k == z \
        and self.check_primitives(basis_i, basis_x):
            return True
        else:
            return False

    def check_primitives(self, basis_i, basis_j):
        primitives_i = basis_i.primitive_gaussian_array
        primitives_j = basis_j.primitive_gaussian_array
        for k in range(len(primitives_i)):
            if primitives_i[k].contraction != primitives_j[k].contraction \
            or primitives_i[k].exponent != primitives_j[k].exponent:
                return False
        return True

    def sort_index(self, a, b, c, d):
        if a > b:
            a, b = b, a
        if c > d:
            c, d = d, c
        if a > c:
            a, c = c, a
            b, d = d, b
        if a == c and b > d:
            a, c = c, a
            b, d = d, b
        return a, b, c, d
