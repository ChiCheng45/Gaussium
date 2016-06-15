from src.main.common import Vector
from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import Molecule
from math import pi
import numpy as np
import copy, heapq


class MoleculeFactory:

    @classmethod
    def point_group(cls, nuclei_array):
        # if only one nucleus then center the nuclei and skip symmetry checks
        if len(nuclei_array) == 1:
            nuclei_array[0].coordinates = (0.0, 0.0, 0.0)
            return nuclei_array

        nuclei_array, rotation, reflection = cls.standard_orientation(nuclei_array)

        # run molecule through the flow diagram
        if cls.check_linear(nuclei_array):

            if cls.check_inversion_symmetry(nuclei_array):
                return Molecule(nuclei_array, rotation, reflection, 'D_{inf h}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'C_{inf v}')

        elif cls.check_high_symmetry(rotation):

            if cls.check_inversion_symmetry(nuclei_array):

                if cls.check_icosahedron(rotation):
                    return Molecule(nuclei_array, rotation, reflection, 'I_{h}')
                else:
                    return Molecule(nuclei_array, rotation, reflection, 'O_{h}')

            else:
                return Molecule(nuclei_array, rotation, reflection, 'T_{d}')

        elif len(rotation) >= 1:
            n = str(cls.get_n_fold(rotation))

            if cls.check_n_two_fold_perpendicular_to_n_fold(rotation):
                pass

            elif cls.check_sigma_h(reflection):
                return Molecule(nuclei_array, rotation, reflection, 'D_{' + n + 'h}')

            elif cls.check_n_sigma_v(n, reflection):
                return Molecule(nuclei_array, rotation, reflection, 'C_{' + n + 'v}')

        elif cls.check_sigma_h(reflection):
            return Molecule(nuclei_array, rotation, reflection, 'C_{s}')

        elif cls.check_inversion_symmetry(nuclei_array):
            return Molecule(nuclei_array, rotation, reflection, 'C_{i}')

        else:
            return Molecule(nuclei_array, rotation, reflection, 'C_{1}')

    @classmethod
    def check_n_sigma_v(cls, n, reflection_symmetry):
        n = int(n)
        i = 0
        for reflection in reflection_symmetry:
            coordinates = Vector.cartesian_to_spherical(reflection.vector)
            if coordinates[1] - pi / 2 <= 1e-3:
                i += 1

        if n == i:
            return True
        else:
            return False

    @classmethod
    def check_n_two_fold_perpendicular_to_n_fold(cls, rotation_symmetry):
        principal_axis = cls.get_principal_axis(rotation_symmetry)

        axis_of_rotation = []
        for rotation in rotation_symmetry:
            if rotation.fold == 2:
                axis_of_rotation.append(rotation.vector)

        vectors = cls.remove_duplicate(cls.cross_products_vertices_vertices(axis_of_rotation)) + [principal_axis.vector]

        spherical_coordinates = []
        for vector in vectors:
            coordinates = Vector.cartesian_to_spherical(vector)
            spherical_coordinates.append(coordinates)

        if all(coordinates[1] % pi <= 1e-3 for coordinates in spherical_coordinates) \
        and all(coordinates[2] <= 1e-3 for coordinates in spherical_coordinates):

            if len(axis_of_rotation) == principal_axis.fold:

                return True

        return False

    @staticmethod
    def get_principal_axis(rotation_symmetry):
        for rotation in rotation_symmetry:
            if Vector.theta(rotation.vector) % pi <= 1e-3:
                return rotation

    @staticmethod
    def check_sigma_h(reflection_symmetry):

        for reflection in reflection_symmetry:
            if Vector.theta(reflection.vector) <= 1e-3:
                return True

        return False

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
            nuclei_array.pop(i)

            if not any(nucleus.element == nuclei.element for nucleus in nuclei_array) \
            or not any(nucleus.coordinates == coordinate_inverse for nucleus in nuclei_array):
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
    def get_n_fold(cls, rotation_symmetry):
        principal_axis = cls.get_principal_axis(rotation_symmetry)
        return principal_axis.fold

    @classmethod
    def standard_orientation(cls, nuclei_array):
        nuclei_array = cls.center_molecule(nuclei_array)
        rotation_symmetry = cls.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry = cls.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry)

        if len(rotation_symmetry) > 1:
            first_highest_symmetry = None
            second_highest_symmetry = None

            highest_n_folds = heapq.nlargest(2, [rotation.fold for rotation in rotation_symmetry])

            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[0]:
                    first_highest_symmetry = rotation
                    break

            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[1] and rotation != first_highest_symmetry:
                    second_highest_symmetry = rotation
                    break

            quaternion = cls.quaternion_rotate_from_phi(first_highest_symmetry.vector, 0.0)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            quaternion = cls.quaternion_to_z_axis(second_highest_symmetry.vector)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 1 and len(reflection_symmetry) >= 1:
            first_highest_symmetry = rotation_symmetry[0]
            reflection_d = None
            sigma_d = False

            quaternion = cls.quaternion_to_z_axis(first_highest_symmetry.vector)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            for reflection in reflection_symmetry:
                if Vector.phi(reflection.vector) > 1e-3:
                    reflection_d = reflection
                    sigma_d = True
                    break

            if sigma_d:
                quaternion = cls.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
                cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 0 and len(reflection_symmetry) > 1:
            reflection_h = reflection_symmetry[0]
            reflection_d = reflection_symmetry[1]

            quaternion = cls.quaternion_to_z_axis(reflection_h.vector)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            quaternion = cls.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 0 and len(reflection_symmetry) == 1:
            quaternion = cls.quaternion_to_z_axis(reflection_symmetry[0].vector)
            cls.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        return nuclei_array, rotation_symmetry, reflection_symmetry

    @staticmethod
    def quaternion_rotate_from_phi(symmetry_vector, angle):
        vector = (0.0, 0.0, 1.0)
        phi = - Vector.phi(symmetry_vector) + angle
        quaternion = Vector.create_quaternion(vector, phi)
        return quaternion

    @staticmethod
    def quaternion_to_z_axis(symmetry_vector):
        vector = (- symmetry_vector[1], symmetry_vector[0], 0.0)
        theta = - Vector.theta(symmetry_vector)
        quaternion = Vector.create_quaternion(vector, theta)
        return quaternion

    @staticmethod
    def rotate_all_vectors(quaternion, rotation_symmetry_list, reflection_symmetry_list, nuclei_array):
        for rotation in rotation_symmetry_list:
            rotation.vector = Vector.quaternion_rotation(quaternion, rotation.vector)

        for reflection in reflection_symmetry_list:
            reflection.vector = Vector.quaternion_rotation(quaternion, reflection.vector)

        for nuclei in nuclei_array:
            nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

    @classmethod
    def brute_force_reflection_symmetry(cls, nuclei_array, rotation_symmetry):
        nuclei_array = cls.remove_center_nuclei(nuclei_array)
        vertices = cls.vertices(nuclei_array)
        edge_center = cls.center_two_vertices(vertices)
        cross_vertices_vertices = cls.cross_products_vertices_vertices(vertices)
        cross_edge_center_vertices = cls.cross_product_center_edge_vertices(edge_center, vertices)

        vector_cross = cls.remove_duplicate(cross_vertices_vertices + cross_edge_center_vertices)

        reflection_planes = []
        if len(vector_cross) > 0:

            # rotate all orthogonal vectors by principal axis'
            vectors_cross_rotated = []
            if len(rotation_symmetry) > 0:

                highest_symmetries = []
                highest_n_fold = heapq.nlargest(1, [rotation.fold for rotation in rotation_symmetry])[0]
                for rotation in rotation_symmetry:
                    if rotation.fold == highest_n_fold:
                        highest_symmetries.append(rotation)

                for highest_symmetry in highest_symmetries:

                    quaternion_list = []
                    vector = highest_symmetry.vector
                    for i in range(1, highest_symmetry.fold):
                        theta = pi * i / highest_symmetry.fold
                        quaternion = Vector.create_quaternion(vector, theta)
                        quaternion_list.append(quaternion)

                    for quaternion in quaternion_list:
                        for orthogonal_vector_i in vector_cross:
                            orthogonal_vector_j = Vector.quaternion_rotation(quaternion, orthogonal_vector_i)
                            vectors_cross_rotated.append(orthogonal_vector_j)

            total_vectors_cross = cross_vertices_vertices + vectors_cross_rotated
            vectors_reflection_plane = cls.remove_duplicate(total_vectors_cross)

            # create householder matrices
            householder_matrices = []
            for planes in vectors_reflection_plane:
                planes = np.matrix(planes)
                householder_matrices.append(np.identity(3) - 2 * planes.T * planes)

            # check each reflection
            for i, matrix in enumerate(householder_matrices):
                if cls.check_reflection(nuclei_array, matrix):
                    planes_of_reflection = ReflectionSymmetry(vectors_reflection_plane[i])
                    reflection_planes.append(planes_of_reflection)

        return reflection_planes

    @classmethod
    def brute_force_rotation_symmetry(cls, nuclei_array):
        nuclei_array = cls.remove_center_nuclei(nuclei_array)
        vertices = cls.vertices(nuclei_array)
        edge_center = cls.center_two_vertices(vertices)
        faces_center = cls.center_three_vertices(vertices)
        cross_vertices_vertices = cls.cross_products_vertices_vertices(vertices)
        cross_edge_center_vertices = cls.cross_product_center_edge_vertices(edge_center, vertices)

        axis_of_rotations_i = cls.remove_duplicate(vertices + edge_center + faces_center + cross_vertices_vertices
        + cross_edge_center_vertices)

        axis_of_rotations_j = []
        if len(axis_of_rotations_i) > 0:

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
                    if cls.check_quaternion(nuclei_array, quaternion_matrix[i, j]):
                        n_fold_symmetry_i[j] = i + 2

            # create the rotation symmetry object and return them if the symmetry if > 1-fold

            for i, symmetry in enumerate(n_fold_symmetry_i):
                if symmetry != 1:
                    rotation_symmetry = RotationSymmetry(n_fold_symmetry_i[i], axis_of_rotations_i[i])
                    axis_of_rotations_j.append(rotation_symmetry)

        return axis_of_rotations_j

    @staticmethod
    def check_reflection(nuclei_array, householder_matrix):
        nuclei_array_copy = copy.deepcopy(nuclei_array)

        for nuclei in nuclei_array_copy:
            coordinates = householder_matrix * np.matrix(nuclei.coordinates).T
            coordinates = coordinates.T.tolist()[0]
            nuclei.coordinates = tuple(coordinates)

        for i, nuclei_i in enumerate(nuclei_array):
            for j, nuclei_j in enumerate(nuclei_array_copy):

                if Vector.distance(nuclei_i.coordinates, nuclei_j.coordinates) <= 1e-3 \
                and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                    break

                if j == len(nuclei_array_copy) - 1:
                    return False

        return True

    @staticmethod
    def check_quaternion(nuclei_array, quaternion):
        nuclei_array_copy = copy.deepcopy(nuclei_array)

        for nuclei in nuclei_array_copy:
            nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)

        for i, nuclei_i in enumerate(nuclei_array):
            for j, nuclei_j in enumerate(nuclei_array_copy):

                if Vector.distance(nuclei_i.coordinates, nuclei_j.coordinates) <= 1e-3 \
                and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                    break

                if j == len(nuclei_array_copy) - 1:
                    return False

        return True

    @staticmethod
    def remove_center_nuclei(nuclei_array):
        for i, nuclei in enumerate(nuclei_array):
            if Vector.rho(nuclei.coordinates) <= 1e-3:
                nuclei_array.pop(i)
        return nuclei_array

    @staticmethod
    def vertices(nuclei_array):
        vertices = []
        for nuclei in nuclei_array:
            coordinates = Vector.normalize(nuclei.coordinates)
            vertices.append(coordinates)
        return vertices

    @staticmethod
    def center_two_vertices(vertices):
        center_of_edge = []
        for axis_i in vertices:
            for axis_j in vertices:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_edge = Vector.add(axis_i, axis_j)
                    if Vector.rho(axis_edge) > 1e-3:
                        axis_edge = Vector.normalize(axis_edge)
                        center_of_edge.append(axis_edge)
        return center_of_edge

    @staticmethod
    def center_three_vertices(vertices):
        center_of_faces = []
        for axis_i in vertices:
            for axis_j in vertices:
                for axis_k in vertices:
                    if Vector.distance(axis_i, axis_j) > 1e-3 and Vector.distance(axis_i, axis_k) > 1e-3 \
                    and Vector.distance(axis_j, axis_k) > 1e-3:
                        axis_face = Vector.add(Vector.add(axis_i, axis_j), axis_k)
                        if Vector.rho(axis_face) > 1e-3:
                            axis_face = Vector.normalize(axis_face)
                            center_of_faces.append(axis_face)
        return center_of_faces

    @staticmethod
    def cross_products_vertices_vertices(vertices):
        cross_products = []
        for axis_i in vertices:
            for axis_j in vertices:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_cross = Vector.cross(axis_i, axis_j)
                    if Vector.rho(axis_cross) > 1e-3:
                        axis_cross = Vector.normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

    @classmethod
    def cross_product_center_edge_vertices(cls, edges, vertices):
        cross_products = []
        for axis_i in edges:
            for axis_j in vertices:
                if Vector.distance(axis_i, axis_j) > 1e-3:
                    axis_cross = Vector.cross(axis_i, axis_j)
                    if Vector.rho(axis_cross) > 1e-3:
                        axis_cross = Vector.normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

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
            if Vector.rho(Vector.add(axis, axis_i)) <= 1e-3 \
            or Vector.rho(Vector.minus(axis, axis_i)) <= 1e-3:
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
