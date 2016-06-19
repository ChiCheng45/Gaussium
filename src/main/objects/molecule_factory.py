from src.main.common import rho, theta, phi
from src.main.common import vector_add
from src.main.common import normalize
from src.main.common import cross_product
from src.main.common import coordinate_distance
from src.main.common import cartesian_to_spherical
from src.main.common import create_quaternion
from src.main.common import quaternion_rotation
from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import Molecule
from src.main.objects import Nuclei
from math import pi
import numpy as np
import copy, heapq


class MoleculeFactory:

    def __init__(self, error=1e-3):
        self.error = error

    def point_group(self, nuclei_array):
        # run molecule through the flow diagram
        nuclei_array = self.center_molecule(nuclei_array)
        if len(nuclei_array) == 1:                              # one nuclei
            return Molecule(nuclei_array, [], [], 'C_{1}')

        rotation, reflection = self.brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = self.standard_orientation(nuclei_array, rotation, reflection)
        if self.check_linear(nuclei_array):                          # Linear
            if self.check_inversion_symmetry(nuclei_array):
                return Molecule(nuclei_array, rotation, reflection, 'D_{inf h}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'C_{inf v}')

        if self.check_high_symmetry(rotation):                       # Polyhedral
            if not self.check_inversion_symmetry(nuclei_array):
                return Molecule(nuclei_array, rotation, reflection, 'T_{d}')
            elif any([vector.fold == 5 for vector in rotation]):
                return Molecule(nuclei_array, rotation, reflection, 'I_{h}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'O_{h}')

        if len(rotation) == 0:                                      # Nonaxial
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, rotation, reflection, 'C_{s}')
            elif self.check_inversion_symmetry(nuclei_array):
                return Molecule(nuclei_array, rotation, reflection, 'C_{i}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'C_{1}')

        n = self.return_principal_axis(rotation).fold
        if self.check_n_two_fold_perpendicular_to_n_fold(rotation):  # Dihedral
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + 'h}')
            elif self.check_n_sigma_v(n, reflection):
                return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + 'd}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + '}')

        else:                                                       # Cyclic
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, rotation, reflection, 'C_{' + str(n) + 'h}')
            elif self.check_n_sigma_v(n, reflection):
                return Molecule(nuclei_array, rotation, reflection, 'C_{' + str(n) + 'v}')
            else:
                return Molecule(nuclei_array, rotation, reflection, 'C_{' + str(n) + '}')

    def check_n_sigma_v(self, n, reflection_symmetry):
        i = 0
        for reflection in reflection_symmetry:
            coordinates = cartesian_to_spherical(reflection.vector)
            if coordinates[1] - pi / 2 <= self.error:
                i += 1
        if n == i:
            return True
        else:
            return False

    def check_n_two_fold_perpendicular_to_n_fold(self, rotation_symmetry):
        principal_axis = self.return_principal_axis(rotation_symmetry)

        axis_of_rotation = []
        for rotation in rotation_symmetry:
            if principal_axis != rotation and rotation.fold == 2:
                axis_of_rotation.append(rotation.vector)

        vectors = self.remove_duplicate(self.cross_products_vertices_vertices(axis_of_rotation)) + [principal_axis.vector]

        spherical_coordinates = []
        for vector in vectors:
            coordinates = cartesian_to_spherical(vector)
            spherical_coordinates.append(coordinates)

        if all(coordinates[1] % pi <= self.error for coordinates in spherical_coordinates) \
        and all(coordinates[2] <= self.error for coordinates in spherical_coordinates):

            if len(axis_of_rotation) == principal_axis.fold:
                return True

        return False

    def check_sigma_h(self, reflection_symmetry):
        for reflection in reflection_symmetry:
            if theta(reflection.vector) % pi <= self.error:
                return True
        return False

    def check_high_symmetry(self, rotation_symmetry):
        i = 0
        for rotation in rotation_symmetry:
            if rotation.fold > 2:
                i += 1
        if i > 2:
            return True
        else:
            return False

    def check_inversion_symmetry(self, nuclei_array):

        nuclei_array_inverse = []
        for nuclei in nuclei_array:
            coordinate_inverse = (-nuclei.coordinates[0], -nuclei.coordinates[1], -nuclei.coordinates[2])
            nuclei_inverse = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinate_inverse)
            nuclei_array_inverse.append(nuclei_inverse)

        for nuclei in nuclei_array:
            for i, nuclei_inverse in enumerate(nuclei_array_inverse):
                if coordinate_distance(nuclei.coordinates, nuclei_inverse.coordinates) <= self.error \
                and (nuclei.charge - nuclei_inverse.charge) == 0.0:
                    break
                if i == len(nuclei_array_inverse) - 1:
                    return False
        return True

    def check_linear(self, nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        spherical_coordinates = []

        for nuclei in nuclei_array_copy:
            nuclei.coordinates = cartesian_to_spherical(nuclei.coordinates)
            spherical_coordinates.append(nuclei)

        if all(nuclei.coordinates[1] % pi <= self.error for nuclei in spherical_coordinates) \
        and all(nuclei.coordinates[2] <= self.error for nuclei in spherical_coordinates):
            return True
        else:
            return False

    def return_principal_axis(self, rotation_symmetry):
        for rotation in rotation_symmetry:
            if theta(rotation.vector) % pi <= self.error:
                return rotation

    def standard_orientation(self, nuclei_array, rotation_symmetry, reflection_symmetry):

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

            quaternion = self.quaternion_rotate_to_z_axis(first_highest_symmetry.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            quaternion = self.quaternion_rotate_from_phi(second_highest_symmetry.vector, 0.0)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 1 and len(reflection_symmetry) >= 1:
            first_highest_symmetry = rotation_symmetry[0]
            reflection_d = None
            sigma_d = False

            quaternion = self.quaternion_rotate_to_z_axis(first_highest_symmetry.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            for reflection in reflection_symmetry:
                if phi(reflection.vector) > self.error:
                    reflection_d = reflection
                    sigma_d = True
                    break

            if sigma_d:
                quaternion = self.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
                self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 0 and len(reflection_symmetry) > 1:
            reflection_h = reflection_symmetry[0]
            reflection_d = reflection_symmetry[1]

            quaternion = self.quaternion_rotate_to_z_axis(reflection_h.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            quaternion = self.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        elif len(rotation_symmetry) == 0 and len(reflection_symmetry) == 1:
            quaternion = self.quaternion_rotate_to_z_axis(reflection_symmetry[0].vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        return nuclei_array, rotation_symmetry, reflection_symmetry

    def rotate_all_vectors(self, quaternion_i, rotation_symmetry_list_i, reflection_symmetry_list_i, nuclei_array_i):
        for rotation_i in rotation_symmetry_list_i:
            rotation_i.vector = quaternion_rotation(quaternion_i, rotation_i.vector)
        for reflection_i in reflection_symmetry_list_i:
            reflection_i.vector = quaternion_rotation(quaternion_i, reflection_i.vector)
        for nuclei_i in nuclei_array_i:
            nuclei_i.coordinates = quaternion_rotation(quaternion_i, nuclei_i.coordinates)

    def quaternion_rotate_from_phi(self, symmetry_vector, angle):
        vector = (0.0, 0.0, 1.0)
        phi_i = - phi(symmetry_vector) + angle
        quaternion = create_quaternion(vector, phi_i)
        return quaternion

    def quaternion_rotate_to_z_axis(self, symmetry_vector):
        vector = (-symmetry_vector[1], symmetry_vector[0], 0.0)
        theta_i = - theta(symmetry_vector)
        quaternion = create_quaternion(vector, theta_i)
        return quaternion

    def brute_force_symmetry(self, nuclei_array):
        nuclei_array = self.remove_center_nuclei(nuclei_array)
        vertices = self.vertices(nuclei_array)
        edge_center = self.center_two_vertices(nuclei_array)
        faces_center = self.center_three_vertices(nuclei_array)
        cross_corners_corners = self.cross_products_vertices_vertices(vertices)
        cross_edge_center_corners = self.cross_product_center_edge_vertices(vertices)

        rotation_symmetry = self.brute_force_rotation_symmetry(nuclei_array, vertices, edge_center, faces_center,
        cross_corners_corners, cross_edge_center_corners)
        reflection_symmetry = self.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry,
        cross_corners_corners, cross_edge_center_corners)

        return rotation_symmetry, reflection_symmetry

    def remove_center_nuclei(self, nuclei_array):
        for i, nuclei in enumerate(nuclei_array):
            if rho(nuclei.coordinates) <= self.error:
                nuclei_array.pop(i)
        return nuclei_array

    def vertices(self, nuclei_array):
        corner = []
        for nuclei in nuclei_array:
            coordinates = normalize(nuclei.coordinates)
            corner.append(coordinates)
        return corner

    def center_two_vertices(self, nuclei_array):
        center_of_edge = []
        for nuclei_i in nuclei_array:
            for nuclei_j in nuclei_array:
                if nuclei_i is not nuclei_j:
                    axis_i = nuclei_i.coordinates
                    axis_j = nuclei_j.coordinates
                    axis_edge = vector_add(axis_i, axis_j)
                    if rho(axis_edge) > self.error:
                        axis_edge = normalize(axis_edge)
                        center_of_edge.append(axis_edge)
        return center_of_edge

    def center_three_vertices(self, nuclei_array):
        center_of_faces = []
        for nuclei_i in nuclei_array:
            for nuclei_j in nuclei_array:
                for nuclei_k in nuclei_array:
                    if (nuclei_i is not nuclei_j) and (nuclei_i is not nuclei_k) and (nuclei_j is not nuclei_k):
                        axis_i = nuclei_i.coordinates
                        axis_j = nuclei_j.coordinates
                        axis_k = nuclei_k.coordinates
                        axis_face = vector_add(vector_add(axis_i, axis_j), axis_k)
                        if rho(axis_face) > self.error:
                            axis_face = normalize(axis_face)
                            center_of_faces.append(axis_face)
        return center_of_faces

    def cross_products_vertices_vertices(self, vertices):
        vertices = self.remove_duplicate(vertices)
        cross_products = []
        for axis_i in vertices:
            for axis_j in vertices:
                if axis_i is not axis_j:
                    axis_cross = cross_product(axis_i, axis_j)
                    if rho(axis_cross) > self.error:
                        axis_cross = normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

    def cross_product_center_edge_vertices(self, vertices):
        vertices = self.remove_duplicate(vertices)
        cross_products = []
        for axis_i in vertices:
            for axis_j in vertices:
                for axis_k in vertices:
                    if (axis_i is not axis_j) and (axis_i is not axis_k) and (axis_j is not axis_k):
                        axis_l = vector_add(axis_i, axis_j)
                        axis_cross = cross_product(axis_l, axis_k)
                        if rho(axis_cross) > self.error:
                            axis_cross = normalize(axis_cross)
                            cross_products.append(axis_cross)
        return cross_products

    def brute_force_reflection_symmetry(self, nuclei_array, rotation_symmetry, cross_corners_corners,
        cross_edge_center_corners):

        # rotate all orthogonal vectors by principal axis by twice it's n-fold
        vector_cross = self.remove_duplicate(cross_corners_corners + cross_edge_center_corners)
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
                    theta_i = pi * i / highest_symmetry.fold
                    quaternion = create_quaternion(vector, theta_i)
                    quaternion_list.append(quaternion)

                for quaternion in quaternion_list:
                    for orthogonal_vector_i in vector_cross:
                        orthogonal_vector_j = quaternion_rotation(quaternion, orthogonal_vector_i)
                        vectors_cross_rotated.append(orthogonal_vector_j)

        reflection_planes = []
        if len(vector_cross) > 0:

            total_vectors_cross = cross_corners_corners + vectors_cross_rotated
            vectors_reflection_plane = self.remove_duplicate(total_vectors_cross)

            # create householder matrices
            householder_matrices = []
            for planes in vectors_reflection_plane:
                planes = np.matrix(planes)
                householder_matrices.append(np.identity(3) - 2 * planes.T * planes)

            # check each reflection
            for i, matrix in enumerate(householder_matrices):
                if self.check_reflection(nuclei_array, matrix):
                    planes_of_reflection = ReflectionSymmetry(vectors_reflection_plane[i])
                    reflection_planes.append(planes_of_reflection)

        return reflection_planes

    def check_reflection(self, nuclei_array, householder_matrix):
        nuclei_array_copy = []
        for nuclei in nuclei_array:
            coordinates = householder_matrix * np.matrix(nuclei.coordinates).T
            coordinates = tuple(coordinates.T.tolist()[0])
            nuclei_copy = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinates)
            nuclei_array_copy.append(nuclei_copy)
        for nuclei_i in nuclei_array:
            for j, nuclei_j in enumerate(nuclei_array_copy):
                if coordinate_distance(nuclei_i.coordinates, nuclei_j.coordinates) <= self.error \
                and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                    break
                if j == len(nuclei_array_copy) - 1:
                    return False
        return True

    def brute_force_rotation_symmetry(self, nuclei_array, corner, edge_center, faces_center, cross_corners_corners,
        cross_edge_corners_vertices):

        axis_of_rotations_i = self.remove_duplicate(corner + edge_center + faces_center + cross_corners_corners
        + cross_edge_corners_vertices)
        axis_of_rotations_j = []
        if len(axis_of_rotations_i) > 0:

            # create quaternion for angle for pi to pi / 4 around the axis of rotation
            quaternion_matrix = np.empty((7, len(axis_of_rotations_i)), dtype=tuple)
            for i in range(7):
                angle = 2 * pi / (i + 2)
                for j, axis in enumerate(axis_of_rotations_i):
                    quaternion_matrix[i, j] = create_quaternion(axis, angle)

            # test all quaternions and create a list of the highest fold symmetry for a given axis
            n_fold_symmetry_i = [1] * len(axis_of_rotations_i)
            for i in range(7):
                for j in range(len(axis_of_rotations_i)):
                    if self.check_rotation(nuclei_array, quaternion_matrix[i, j]):
                        n_fold_symmetry_i[j] = i + 2

            # create the rotation symmetry object and return them if the symmetry if > 1-fold
            for i, symmetry in enumerate(n_fold_symmetry_i):
                if symmetry != 1:
                    rotation_symmetry = RotationSymmetry(n_fold_symmetry_i[i], axis_of_rotations_i[i])
                    axis_of_rotations_j.append(rotation_symmetry)

        return axis_of_rotations_j

    def check_rotation(self, nuclei_array, quaternion):
        nuclei_array_copy = []
        for nuclei in nuclei_array:
            coordinates = quaternion_rotation(quaternion, nuclei.coordinates)
            nuclei_copy = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinates)
            nuclei_array_copy.append(nuclei_copy)
        for nuclei_i in nuclei_array:
            for k, nuclei_k in enumerate(nuclei_array_copy):
                if coordinate_distance(nuclei_i.coordinates, nuclei_k.coordinates) <= self.error \
                and (nuclei_i.charge - nuclei_k.charge) == 0.0:
                    break
                if k == len(nuclei_array_copy) - 1:
                    return False
        return True

    def remove_duplicate(self, axis_of_rotations):
        axis_of_rotations_i = []
        for axis_i in axis_of_rotations:
            if self.check_array(axis_i, axis_of_rotations_i):
                axis_of_rotations_i.append(axis_i)
        return axis_of_rotations_i

    def check_array(self, axis, axis_of_rotations_i):
        for axis_i in axis_of_rotations_i:
            if rho(cross_product(axis, axis_i)) <= self.error:
                return False
        return True

    def center_molecule(self, nuclei_array):
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

        nuclei_array.sort(key=lambda nucleus: rho(nucleus.coordinates))

        if rho(nuclei_array[0].coordinates) <= self.error:

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
