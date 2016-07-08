from src.main.common import rho, theta, phi
from src.main.common import vector_add
from src.main.common import normalize
from src.main.common import cross_product
from src.main.common import coordinate_distance
from src.main.common import cartesian_to_spherical
from src.main.common import create_quaternion
from src.main.common import quaternion_rotation
from src.main.common import create_householder_matrix
from src.main.common import householder_matrix_reflection
from src.main.objects import PointGroup
from src.main.objects import D4h
from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import InversionSymmetry
from src.main.objects import Molecule
from src.main.objects import Nuclei
from math import pi
import numpy as np
import copy, heapq


class MoleculeFactory:

    def __init__(self, error=1e-2):
        self.error = error

    def point_group(self, nuclei_array):
        nuclei_array = self.center_molecule(nuclei_array)

        if len(nuclei_array) == 1:                                      # Point
            return Molecule(nuclei_array, D4h())

        rotation, reflection, inversion = self.brute_force_symmetry(nuclei_array)
        self.standard_orientation(nuclei_array, rotation, reflection)

        if self.check_linear(nuclei_array):                             # Linear
            rotation = reflection = []
            if len(inversion) == 1:
                return Molecule(nuclei_array, D4h())
            else:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{inf v}'))

        if self.check_high_symmetry(rotation):                           # Polyhedral
            if not len(inversion) == 1:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'T_{d}'))
            elif any([vector.fold == 5 for vector in rotation]):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'I_{h}'))
            else:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'O_{h}'))

        if len(rotation) == 0:                                          # Nonaxial
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{s}'))
            elif len(inversion) == 1:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{i}'))
            else:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{1}'))

        n = self.return_principal_axis(rotation).fold

        if self.check_n_two_fold_perpendicular_to_n_fold(rotation):     # Dihedral
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'D_{' + str(n) + 'h}'))
            elif self.check_n_sigma_v(n, reflection):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'D_{' + str(n) + 'd}'))
            else:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'D_{' + str(n) + '}'))

        else:                                                           # Cyclic
            if self.check_sigma_h(reflection):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{' + str(n) + 'h}'))
            elif self.check_n_sigma_v(n, reflection):
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{' + str(n) + 'v}'))
            else:
                return Molecule(nuclei_array, PointGroup(rotation, reflection, inversion, 'C_{' + str(n) + '}'))

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
            if principal_axis != rotation and rotation.fold == 2 and (theta(rotation.vector) - pi/2 <= self.error):
                axis_of_rotation.append(rotation.vector)

        if len(axis_of_rotation) == principal_axis.fold:
            return True
        else:
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

    def check_linear(self, nuclei_array):
        for nuclei in nuclei_array:
            if theta(nuclei.coordinates) % pi > self.error:
                return False
        return True

    def return_principal_axis(self, rotation_symmetry):
        for rotation in rotation_symmetry:
            if theta(rotation.vector) % pi <= self.error:
                return rotation

    def standard_orientation(self, nuclei_array, rotation_symmetry, reflection_symmetry):

        if len(rotation_symmetry) > 1:
            highest_n_folds = heapq.nlargest(2, [rotation.fold for rotation in rotation_symmetry])

            first_highest_symmetry = None
            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[0]:
                    first_highest_symmetry = rotation
                    break

            second_highest_symmetry = None
            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[1] and rotation != first_highest_symmetry:
                    second_highest_symmetry = rotation
                    break

            quaternion = self.quaternion_rotate_to_z_axis(first_highest_symmetry.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)
            quaternion = self.quaternion_rotate_from_phi(second_highest_symmetry.vector, 0.0)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        if len(rotation_symmetry) == 1 and len(reflection_symmetry) >= 1:
            first_highest_symmetry = rotation_symmetry[0]

            quaternion = self.quaternion_rotate_to_z_axis(first_highest_symmetry.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

            reflection_d = None
            for reflection in reflection_symmetry:
                if phi(reflection.vector) > self.error:
                    reflection_d = reflection
                    break

            if reflection_d is not None:
                quaternion = self.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
                self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) > 1:
            reflection_h = reflection_symmetry[0]
            reflection_d = reflection_symmetry[1]

            quaternion = self.quaternion_rotate_to_z_axis(reflection_h.vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)
            quaternion = self.quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        if rho(nuclei_array[0].coordinates) > 1e-3:
            i = 0
            j = 1
        else:
            i = 1
            j = 2

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) == 1:
            quaternion = self.quaternion_rotate_to_z_axis(reflection_symmetry[0].vector)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)
            quaternion = self.quaternion_rotate_from_phi(nuclei_array[i].coordinates, 0.0)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) == 0:
            quaternion = self.quaternion_rotate_to_z_axis(nuclei_array[i].coordinates)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)
            quaternion = self.quaternion_rotate_from_phi(nuclei_array[j].coordinates, 0.0)
            self.rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

    def rotate_all_vectors(self, quaternion, rotation_symmetry_list, reflection_symmetry_list, nuclei_array):
        for rotation in rotation_symmetry_list:
            rotation.vector = quaternion_rotation(quaternion, rotation.vector)
        for reflection in reflection_symmetry_list:
            reflection.vector = quaternion_rotation(quaternion, reflection.vector)
        for nuclei in nuclei_array:
            nuclei.coordinates = quaternion_rotation(quaternion, nuclei.coordinates)

    def quaternion_rotate_from_phi(self, symmetry_vector, angle):
        vector = (0.0, 0.0, 1.0)
        phi_i = - phi(symmetry_vector) + angle
        quaternion = create_quaternion(vector, phi_i)
        return quaternion

    def quaternion_rotate_to_z_axis(self, symmetry_vector):
        vector = (-symmetry_vector[1], symmetry_vector[0], 0.0)
        if rho(vector) <= self.error:
            vector = (1.0, 0.0, 0.0)
        theta_i = - theta(symmetry_vector)
        quaternion = create_quaternion(vector, theta_i)
        return quaternion

    def brute_force_symmetry(self, nuclei_array):
        nuclei_array = self.remove_center_nuclei(nuclei_array)
        vertices = self.remove_duplicate(self.vertices(nuclei_array))
        edge_center = self.remove_duplicate(self.center_two_vertices(nuclei_array))
        cross_vertices_vertices = self.cross_products_vertices_vertices(vertices)
        cross_edge_vertices = self.cross_product_edge_vertices(vertices, edge_center)
        cross_edge_edge = self.cross_product_edge_edge(edge_center)

        rotation_symmetry = self.brute_force_rotation_symmetry(nuclei_array, vertices, edge_center,
        cross_vertices_vertices, cross_edge_vertices, cross_edge_edge)

        reflection_symmetry = self.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry,
        vertices, cross_vertices_vertices, cross_edge_vertices)

        inversion_symmetry = InversionSymmetry()
        if self.check_symmetry_operation(nuclei_array, inversion_symmetry):
            inversion_symmetry = [inversion_symmetry]
        else:
            inversion_symmetry = []

        return rotation_symmetry, reflection_symmetry, inversion_symmetry

    def remove_center_nuclei(self, nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        for i, nuclei in enumerate(nuclei_array_copy):
            if rho(nuclei.coordinates) <= self.error:
                nuclei_array_copy.pop(i)
        return nuclei_array_copy

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

    def cross_products_vertices_vertices(self, vertices):
        cross_products = []
        for axis_i in vertices:
            for axis_j in vertices:
                if axis_i is not axis_j:
                    axis_cross = cross_product(axis_i, axis_j)
                    if rho(axis_cross) > self.error:
                        axis_cross = normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

    def cross_product_edge_vertices(self, vertices, center_of_edge):
        cross_products = []
        for axis_i in vertices:
            for axis_j in center_of_edge:
                    axis_cross = cross_product(axis_i, axis_j)
                    if rho(axis_cross) > self.error:
                        axis_cross = normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

    def cross_product_edge_edge(self, center_of_edge):
        cross_products = []
        for axis_i in center_of_edge:
            for axis_j in center_of_edge:
                if axis_i is not axis_j:
                    axis_cross = cross_product(axis_i, axis_j)
                    if rho(axis_cross) > self.error:
                        axis_cross = normalize(axis_cross)
                        cross_products.append(axis_cross)
        return cross_products

    def brute_force_reflection_symmetry(self, nuclei_array, rotation_symmetry, vertices, cross_vertices_vertices,
        cross_edge_vertices):

        # rotate all orthogonal vectors by principal axis by twice it's n-fold
        vector_cross = self.remove_duplicate(vertices + cross_vertices_vertices + cross_edge_vertices)
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

            total_vectors_cross = cross_vertices_vertices + vectors_cross_rotated
            vectors_reflection_plane = self.remove_duplicate(total_vectors_cross)

            planes_of_reflection = []
            for planes in vectors_reflection_plane:
                householder_matrix = ReflectionSymmetry(planes)
                planes_of_reflection.append(householder_matrix)

            # check each reflection
            for plane_of_reflection in planes_of_reflection:
                if self.check_symmetry_operation(nuclei_array, plane_of_reflection):
                    reflection_planes.append(plane_of_reflection)

        return reflection_planes

    def brute_force_rotation_symmetry(self, nuclei_array, corner, edge_center, cross_vertices_vertices,
        cross_edge_vertices, cross_edge_edge):

        axis_of_rotations_i = self.remove_duplicate(corner + edge_center + cross_vertices_vertices
        + cross_edge_vertices + cross_edge_edge)
        axis_of_rotations_j = []
        if len(axis_of_rotations_i) > 0:

            # create quaternion for angle for pi to pi / 4 around the axis of rotation
            quaternion_matrix = np.empty((7, len(axis_of_rotations_i)), dtype=object)
            for i in range(7):
                for j, axis in enumerate(axis_of_rotations_i):
                    quaternion_matrix[i, j] = RotationSymmetry(i + 2, axis)

            # test all quaternions and create a list of the highest fold symmetry for a given axis
            n_fold_symmetry_i = [1] * len(axis_of_rotations_i)
            for i in range(7):
                for j in range(len(axis_of_rotations_i)):
                    if self.check_symmetry_operation(nuclei_array, quaternion_matrix[i, j]):
                        n_fold_symmetry_i[j] = i + 2

            # create the rotation symmetry object and return them if the symmetry if > 1-fold
            for i, symmetry in enumerate(n_fold_symmetry_i):
                if symmetry != 1:
                    rotation_symmetry = quaternion_matrix[n_fold_symmetry_i[i] - 2, i]
                    axis_of_rotations_j.append(rotation_symmetry)

        return axis_of_rotations_j

    def check_symmetry_operation(self, nuclei_array, symmetry):
        nuclei_array_copy = []
        for nuclei in nuclei_array:
            coordinates = symmetry.operate(nuclei.coordinates)
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
