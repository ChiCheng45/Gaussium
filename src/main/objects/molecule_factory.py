from src.main.common import Vector
from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import Molecule
from math import pi
import numpy as np
import copy, heapq


class MoleculeFactory:

    def create(self, nuclei_array):
        pass

    @classmethod
    def point_group(cls, nuclei_array):
        number_of_nuclei = len(nuclei_array)

        # if only one nucleus then center the nuclei and skip symmetry checks
        if number_of_nuclei == 1:
            nuclei_array[0].coordinates = (0.0, 0.0, 0.0)
            return nuclei_array

        else:

            # put the molecule in the standard orientation
            nuclei_array, rotation_symmetry = cls.standard_orientation(nuclei_array)

            # run molecule through the flow diagram
            if cls.check_linear(nuclei_array):

                if cls.check_inversion_symmetry(nuclei_array):

                    return Molecule(nuclei_array, 'D_{inf h}')
                else:
                    return Molecule(nuclei_array, 'C_{inf v}')

            elif cls.check_high_symmetry(rotation_symmetry):

                if cls.check_inversion_symmetry(nuclei_array):

                    if cls.check_icosahedron(rotation_symmetry):
                        return Molecule(nuclei_array, 'I_{h}')
                    else:
                        return Molecule(nuclei_array, 'O_{h}')

                else:
                    return Molecule(nuclei_array, 'T_{d}')

            elif len(rotation_symmetry) >= 1:
                pass

            else:
                pass

    @staticmethod
    def check_icosahedron(rotation_symmetry):

        if any([vector.fold == 5 for vector in rotation_symmetry]):
            return True
        else:
            return False

    @staticmethod
    def check_high_symmetry(rotation_symmetry):
        i = 0
        for rotation in rotation_symmetry:
            if rotation.fold > 2:
                i += 1
        if i > 2:
            return True
        else:
            return False

    @staticmethod
    def check_inversion_symmetry(nuclei_array):

        for i, nuclei in enumerate(nuclei_array):
            coordinate_inverse = (- nuclei.coordinates[0], - nuclei.coordinates[1], - nuclei.coordinates[2])
            check_list = copy.deepcopy(nuclei_array)
            check_list.pop(i)

            if not any(nucleus.element == nuclei.element for nucleus in check_list) \
            or not any(nucleus.coordinates == coordinate_inverse for nucleus in check_list):

                return False

        return True

    @staticmethod
    def check_linear(nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        spherical_coordinates = []

        for nuclei in nuclei_array_copy:
            nuclei.coordinates = Vector.cartesian_to_spherical(nuclei.coordinates)
            spherical_coordinates.append(nuclei)

        if all(nuclei.coordinates[1] % pi <= 1e-3 for nuclei in spherical_coordinates) \
        and all(nuclei.coordinates[2] <= 1e-3 for nuclei in spherical_coordinates):
            return True
        else:
            return False

    @classmethod
    def standard_orientation(cls, nuclei_array):
        nuclei_array = cls.center_molecule(nuclei_array)
        rotation_symmetry_list = cls.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry_list = cls.brute_force_reflection_symmetry(nuclei_array)

        if len(rotation_symmetry_list) > 1:
            first_highest_symmetry = None
            second_highest_symmetry = None

            highest_n_folds = heapq.nlargest(2, [rotation.fold for rotation in rotation_symmetry_list])

            for rotation in rotation_symmetry_list:
                if rotation.fold == highest_n_folds[0]:
                    first_highest_symmetry = rotation
                    break

            for rotation in rotation_symmetry_list:
                if rotation.fold == highest_n_folds[1] and rotation != first_highest_symmetry:
                    second_highest_symmetry = rotation
                    break

            vector = (-first_highest_symmetry.vector[1], first_highest_symmetry.vector[0], 0.0)
            theta = -Vector.theta(first_highest_symmetry.vector)
            quaternion = Vector.create_quaternion(vector, theta)

            for rotation in rotation_symmetry_list:
                rotation.vector = Vector.quaternion_rotation(quaternion, rotation.vector)

            for nuclei in nuclei_array:
                nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

            vector = (0.0, 0.0, 1.0)
            theta = -Vector.phi(second_highest_symmetry.vector)
            quaternion = Vector.create_quaternion(vector, theta)

            for rotation in rotation_symmetry_list:
                rotation.vector = Vector.quaternion_rotation(quaternion, rotation.vector)

            for nuclei in nuclei_array:
                nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

            return nuclei_array, rotation_symmetry_list

        elif len(rotation_symmetry_list) == 1:
            first_highest_symmetry = rotation_symmetry_list[0]

            vector = (- first_highest_symmetry.vector[1], first_highest_symmetry.vector[0], 0.0)
            theta = - Vector.theta(first_highest_symmetry.vector)
            quaternion = Vector.create_quaternion(vector, theta)

            for rotation in rotation_symmetry_list:
                rotation.vector = Vector.quaternion_rotation(quaternion, rotation.vector)

            for nuclei in nuclei_array:
                nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

            return nuclei_array, rotation_symmetry_list

        else:
            return nuclei_array, rotation_symmetry_list

    @classmethod
    def brute_force_reflection_symmetry(cls, nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        vectors = []
        vectors_cross = []
        householder_matrices = []

        # remove the nuclei from the center
        for i, nuclei in enumerate(nuclei_array_copy):
            if Vector.rho(nuclei.coordinates) <= 1e-3:
                nuclei_array_copy.pop(i)

        for nuclei in nuclei_array_copy:
            coordinates = Vector.normalize(nuclei.coordinates)
            vectors.append(coordinates)

        # brute force vector perpendicular to two vertices
        for axis_i in vectors:
            for axis_j in vectors:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_cross = Vector.cross(axis_i, axis_j)
                    if Vector.rho(axis_cross) > 1e-3:
                        axis_cross = Vector.normalize(axis_cross)
                        vectors_cross.append(axis_cross)

        vectors_cross = cls.remove_duplicate(vectors_cross)

        # create householder matrices
        for planes in vectors_cross:
            planes = np.matrix(planes)
            householder_matrices.append(np.identity(3) - 2 * planes.T * planes)

        # check each reflection
        reflection_planes = []
        for i, matrix in enumerate(householder_matrices):
            if cls.check_reflection(nuclei_array_copy, matrix):
                planes_of_reflection = ReflectionSymmetry(vectors_cross[i])
                reflection_planes.append(planes_of_reflection)

        print(reflection_planes)

        return []

    @staticmethod
    def check_reflection(nuclei_array, householder_matrix):
        nuclei_array_copy = copy.deepcopy(nuclei_array)

        for nuclei in nuclei_array_copy:
            coordinates = householder_matrix * np.matrix(nuclei.coordinates).T
            coordinates = coordinates.T.tolist()[0]
            nuclei.coordinates = tuple(coordinates)

        for i, nuclei_i in enumerate(nuclei_array):
            for j, nuclei_j in enumerate(nuclei_array_copy):

                if Vector.distance(nuclei_i.coordinates, nuclei_j.coordinates) < 1e-3 \
                and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                    break

                if j == len(nuclei_array_copy) - 1:
                    return False

        return True

    @classmethod
    def brute_force_rotation_symmetry(cls, nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        axis_of_rotations_vertices = []
        axis_of_rotations_edges = []
        axis_of_rotations_faces = []
        axis_of_rotations_cross = []

        # remove the nuclei from the center
        for i, nuclei in enumerate(nuclei_array_copy):
            if Vector.rho(nuclei.coordinates) <= 1e-3:
                nuclei_array_copy.pop(i)

        # add all vertices
        for nuclei in nuclei_array_copy:
            coordinates = Vector.normalize(nuclei.coordinates)
            axis_of_rotations_vertices.append(coordinates)

        # brute force all vectors in-between points
        for axis_i in axis_of_rotations_vertices:
            for axis_j in axis_of_rotations_vertices:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_edge = Vector.add(axis_i, axis_j)
                    if Vector.rho(axis_edge) > 1e-3:
                        axis_edge = Vector.normalize(axis_edge)
                        axis_of_rotations_edges.append(axis_edge)

        # brute force vectors in the middle of faces with odd numbers of edges
        for axis_i in axis_of_rotations_vertices:
            for axis_j in axis_of_rotations_vertices:
                for axis_k in axis_of_rotations_vertices:
                    if Vector.distance(axis_i, axis_j) > 1e-3 and Vector.distance(axis_i, axis_k) > 1e-3 \
                    and Vector.distance(axis_j, axis_k) > 1e-3:
                        axis_face = Vector.add(Vector.add(axis_i, axis_j), axis_k)
                        if Vector.rho(axis_face) > 1e-3:
                            axis_face = Vector.normalize(axis_face)
                            axis_of_rotations_faces.append(axis_face)

        # brute force vector perpendicular to two vertices
        for axis_i in axis_of_rotations_vertices:
            for axis_j in axis_of_rotations_vertices:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_cross = Vector.cross(axis_i, axis_j)
                    if Vector.rho(axis_cross) > 1e-3:
                        axis_cross = Vector.normalize(axis_cross)
                        axis_of_rotations_cross.append(axis_cross)

        axis_of_rotations_i = axis_of_rotations_vertices + axis_of_rotations_edges + axis_of_rotations_faces \
        + axis_of_rotations_cross
        axis_of_rotations_i = cls.remove_duplicate(axis_of_rotations_i)

        # create quaternion for angle for pi to pi / 4 around the axis of rotation
        quaternion_matrix = np.empty((7, len(axis_of_rotations_i)), dtype=tuple)
        for i in range(7):
            angle = 2 * pi / (i + 2)
            for j, axis in enumerate(axis_of_rotations_i):
                quaternion_matrix[i, j] = Vector.create_quaternion(axis, angle)

        # test all quaternions and create a list of the highest fold symmetry for a given axis
        n_fold_symmetry_i = [1] * len(axis_of_rotations_i)
        for i in range(7):
            for j in range(len(axis_of_rotations_i)):
                if cls.check_quaternion(nuclei_array_copy, quaternion_matrix[i, j]):
                    n_fold_symmetry_i[j] = i + 2

        # create the rotation symmetry object and return them if the symmetry if > 1-fold
        axis_of_rotations = []
        for i, symmetry in enumerate(n_fold_symmetry_i):
            if symmetry != 1:
                rotation_symmetry = RotationSymmetry(n_fold_symmetry_i[i], axis_of_rotations_i[i])
                axis_of_rotations.append(rotation_symmetry)

        return axis_of_rotations

    @staticmethod
    def check_quaternion(nuclei_array, quaternion):
        nuclei_array_copy = copy.deepcopy(nuclei_array)

        for nuclei in nuclei_array_copy:
            nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

        for i, nuclei_i in enumerate(nuclei_array):
            for j, nuclei_j in enumerate(nuclei_array_copy):

                if Vector.distance(nuclei_i.coordinates, nuclei_j.coordinates) < 1e-3 \
                and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                    break

                if j == len(nuclei_array_copy) - 1:
                    return False

        return True

    @classmethod
    def remove_duplicate(cls, axis_of_rotations):
        axis_of_rotations_i = []
        for i, axis_i in enumerate(axis_of_rotations):
            if cls.check_array(axis_i, axis_of_rotations_i):
                axis_of_rotations_i.append(axis_i)
        return axis_of_rotations_i

    @staticmethod
    def check_array(axis, axis_of_rotations_i):
        for axis_i in axis_of_rotations_i:
            if Vector.rho(Vector.add(axis, axis_i)) < 1e-3 \
            or Vector.rho(Vector.minus(axis, axis_i)) < 1e-3:
                return False
        return True

    @staticmethod
    def center_molecule(nuclei_array):
        number_of_nuclei = len(nuclei_array)
        center = [0, 0, 0]

        for nuclei in nuclei_array:
            center[0] += nuclei.coordinates[0] / number_of_nuclei
            center[1] += nuclei.coordinates[1] / number_of_nuclei
            center[2] += nuclei.coordinates[2] / number_of_nuclei

        for nuclei in nuclei_array:
            x = nuclei.coordinates[0] - center[0]
            y = nuclei.coordinates[1] - center[1]
            z = nuclei.coordinates[2] - center[2]
            nuclei.coordinates = (x, y, z)

        nuclei_array.sort(key=lambda nucleus: Vector.rho(nucleus.coordinates))

        if Vector.rho(nuclei_array[0].coordinates) <= 1e-3:

            for i, nuclei in enumerate(nuclei_array):
                if i == 0:
                    translation = nuclei.coordinates
                    nuclei.coordinates = (0.0, 0.0, 0.0)
                else:
                    x = nuclei.coordinates[0] - translation[0]
                    y = nuclei.coordinates[1] - translation[1]
                    z = nuclei.coordinates[2] - translation[2]
                    nuclei.coordinates = (x, y, z)

            return nuclei_array

        elif not Vector.rho(nuclei_array[0].coordinates) <= 1e-3:

            return nuclei_array
