import copy, itertools
from math import pi
import numpy as np
from src.common import coordinate_distance
from src.common import create_quaternion
from src.common import cross_product
from src.common import normalize
from src.common import quaternion_rotation
from src.common import rho
from src.common import vector_add
from src.objects import ImproperRotationSymmetry
from src.objects import InversionSymmetry
from src.objects import Nuclei
from src.objects import ReflectionSymmetry
from src.objects import RotationSymmetry


class SymmetryFactory:

    def __init__(self, error=1e-2):
        self.error = error

    def brute_force_symmetry(self, nuclei_array):
        nuclei_array = self.remove_center_nuclei(nuclei_array)

        vertices = self.remove_duplicate(self.vertices(nuclei_array))
        edge_center = self.remove_duplicate(self.center_two_vertices(nuclei_array))
        cross_vertices_vertices = self.cross_products(vertices, vertices)
        cross_edge_vertices = self.cross_products(vertices, edge_center)
        cross_edge_edge = self.cross_products(edge_center, edge_center)

        rotation_symmetry = self.brute_force_rotation_symmetry(
            nuclei_array, vertices, edge_center, cross_vertices_vertices, cross_edge_vertices, cross_edge_edge
        )
        reflection_symmetry = self.brute_force_reflection_symmetry(
            nuclei_array, rotation_symmetry, vertices, cross_vertices_vertices, cross_edge_vertices
        )
        improper_rotation = self.brute_force_improper_rotation(nuclei_array, rotation_symmetry)
        inversion_symmetry = self.brute_force_inversion_symmetry(nuclei_array)

        return rotation_symmetry, reflection_symmetry, improper_rotation, inversion_symmetry

    def brute_force_rotation_symmetry(self, nuclei_array, corner, edge_center, cross_vertices_vertices,
    cross_edge_vertices, cross_edge_edge):

        axis_of_rotations_i = self.remove_duplicate(
            corner + edge_center + cross_vertices_vertices + cross_edge_vertices + cross_edge_edge
        )

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

    def brute_force_reflection_symmetry(self, nuclei_array, rotation_symmetry, vertices, cross_vertices_vertices,
    cross_edge_vertices):

        # rotate all orthogonal vectors by principal axis by half it's n-fold angle
        vector_cross = self.remove_duplicate(vertices + cross_vertices_vertices + cross_edge_vertices)

        vectors_cross_rotated = []
        if len(rotation_symmetry) > 0:

            highest_symmetries = []
            highest_n_fold = max([rotation.fold for rotation in rotation_symmetry])
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

    def brute_force_improper_rotation(self, nuclei_array, rotation_symmetry):
        improper_rotation_list_i = []
        for rotation in rotation_symmetry:
            improper_rotation = ImproperRotationSymmetry(2 * rotation.fold, rotation.vector)
            improper_rotation_list_i.append(improper_rotation)

        improper_rotation_list_j = []
        for improper_rotation in improper_rotation_list_i:
            if self.check_symmetry_operation(nuclei_array, improper_rotation):
                improper_rotation_list_j.append(improper_rotation)
                continue
            improper_rotation.fold //= 2
            if improper_rotation.fold > 2 and self.check_symmetry_operation(nuclei_array, improper_rotation):
                improper_rotation_list_j.append(improper_rotation)

        return improper_rotation_list_j

    def brute_force_inversion_symmetry(self, nuclei_array):
        inversion = InversionSymmetry()

        inversion_symmetry = []
        if self.check_symmetry_operation(nuclei_array, inversion):
            inversion_symmetry.append(inversion)

        return inversion_symmetry

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
        for nuclei_i, nuclei_j in itertools.combinations(nuclei_array, 2):
            axis_i = nuclei_i.coordinates
            axis_j = nuclei_j.coordinates
            axis_edge = vector_add(axis_i, axis_j)
            if rho(axis_edge) > self.error:
                axis_edge = normalize(axis_edge)
                center_of_edge.append(axis_edge)
        return center_of_edge

    def cross_products(self, vector_i, vector_j):
        cross_products = []
        for axis_i, axis_j in itertools.product(vector_i, vector_j):
            if axis_i is not axis_j:
                axis_cross = cross_product(axis_i, axis_j)
                if rho(axis_cross) > self.error:
                    axis_cross = normalize(axis_cross)
                    cross_products.append(axis_cross)
        return cross_products

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
