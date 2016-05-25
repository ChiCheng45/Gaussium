import copy
import numpy as np
from math import pi
from src.main.common.vector_manipulation import Vector
from src.main.objects import Molecule


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
            nuclei_array = cls.standard_orientation(nuclei_array)

            return nuclei_array

    @classmethod
    def standard_orientation(cls, nuclei_array):

        nuclei_array = cls.center_molecule(nuclei_array)
        cls.brute_force_rotation_symmetry(nuclei_array)

        return nuclei_array

    @classmethod
    def brute_force_rotation_symmetry(cls, nuclei_array):
        nuclei_array_copy = copy.deepcopy(nuclei_array)
        axis_of_rotations_vertices = []
        axis_of_rotations_edges = []
        axis_of_rotations_faces = []

        # remove the molecule from the center
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
                    axis_edge = Vector.normalize(axis_edge)
                    if Vector.rho(axis_edge) > 1e-3:
                        axis_of_rotations_edges.append(axis_edge)

        # brute force vectors in the middle of faces with odd numbers of edges
        for axis_i in axis_of_rotations_vertices:
            for axis_j in axis_of_rotations_vertices:
                for axis_k in axis_of_rotations_vertices:
                    if Vector.distance(axis_i, axis_j) > 1e-3 and Vector.distance(axis_i, axis_k) > 1e-3 \
                    and Vector.distance(axis_j, axis_k) > 1e-3:
                        axis_face = Vector.add(Vector.add(axis_i, axis_j), axis_k)
                        axis_face = Vector.normalize(axis_face)
                        if Vector.rho(axis_face) > 1e-3:
                            axis_of_rotations_faces.append(axis_face)

        axis_of_rotations_total = axis_of_rotations_vertices + axis_of_rotations_edges + axis_of_rotations_faces
        axis_of_rotations_total = cls.remove_duplicate(axis_of_rotations_total)

        quaternion_matrix = np.empty((5, len(axis_of_rotations_total)), dtype=tuple)
        for i in range(5):
            angle = 2 * pi / (i + 2)
            for j, axis in enumerate(axis_of_rotations_total):
                quaternion_matrix[i, j] = Vector.create_quaternion(axis, angle)

        n_fold_symmetry = np.zeros(len(axis_of_rotations_total))
        for i in range(5):
            for j in range(len(axis_of_rotations_total)):
                if cls.check_quaternion(nuclei_array_copy, quaternion_matrix[i, j]):
                    n_fold_symmetry[j] = i + 2

        print(n_fold_symmetry)
        print()

        return axis_of_rotations_total

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
    def remove_duplicate(cls, axis_of_rotations_i):
        axis_of_rotations = []
        for i, axis_i in enumerate(axis_of_rotations_i):
            if cls.check_array(axis_i, axis_of_rotations):
                axis_of_rotations.append(axis_i)
        return axis_of_rotations

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
