from src.main.common import Vector
from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import Molecule
from src.main.objects import Nuclei
from math import pi
import numpy as np
import copy, heapq


def point_group(nuclei_array):
    # run molecule through the flow diagram
    nuclei_array = center_molecule(nuclei_array)
    if len(nuclei_array) == 1:                              # one nuclei
        return Molecule(nuclei_array, [], [], 'C_{1}')

    rotation, reflection = brute_force_symmetry(nuclei_array)
    nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
    if check_linear(nuclei_array):                          # Linear
        if check_inversion_symmetry(nuclei_array):
            return Molecule(nuclei_array, rotation, reflection, 'D_{inf h}')
        else:
            return Molecule(nuclei_array, rotation, reflection, 'C_{inf v}')

    if check_high_symmetry(rotation):                       # Polyhedral
        if not check_inversion_symmetry(nuclei_array):
            return Molecule(nuclei_array, rotation, reflection, 'T_{d}')
        elif check_icosahedron(rotation):
            return Molecule(nuclei_array, rotation, reflection, 'I_{h}')
        else:
            return Molecule(nuclei_array, rotation, reflection, 'O_{h}')

    if len(rotation) == 0:                                      # Nonaxial
        if check_sigma_h(reflection):
            return Molecule(nuclei_array, rotation, reflection, 'C_{s}')
        elif check_inversion_symmetry(nuclei_array):
            return Molecule(nuclei_array, rotation, reflection, 'C_{i}')
        else:
            return Molecule(nuclei_array, rotation, reflection, 'C_{1}')

    n = return_principal_axis(rotation).fold
    if check_n_two_fold_perpendicular_to_n_fold(rotation):  # Dihedral
        if check_sigma_h(reflection):
            return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + 'h}')
        elif check_n_sigma_v(n, reflection):
            return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + 'v}')
        else:
            return Molecule(nuclei_array, rotation, reflection, 'D_{' + str(n) + '}')

    else:                                                       # Cyclic
        if check_sigma_h(reflection):
            return Molecule(nuclei_array, rotation, reflection, 'C_{' + str(n) + 'h}')
        elif check_n_sigma_v(n, reflection):
            return Molecule(nuclei_array, rotation, reflection, 'C_{' + str(n) + 'v}')


def check_n_sigma_v(n, reflection_symmetry):
    i = 0
    for reflection in reflection_symmetry:
        coordinates = Vector.cartesian_to_spherical(reflection.vector)
        if coordinates[1] - pi / 2 <= 1e-3:
            i += 1
    if n == i:
        return True
    else:
        return False


def check_n_two_fold_perpendicular_to_n_fold(rotation_symmetry):
    principal_axis = return_principal_axis(rotation_symmetry)

    axis_of_rotation = []
    for rotation in rotation_symmetry:
        if principal_axis != rotation and rotation.fold == 2:
            axis_of_rotation.append(rotation.vector)

    vectors = remove_duplicate(cross_products_vertices_vertices(axis_of_rotation)) + [principal_axis.vector]

    spherical_coordinates = []
    for vector in vectors:
        coordinates = Vector.cartesian_to_spherical(vector)
        spherical_coordinates.append(coordinates)

    if all(coordinates[1] % pi <= 1e-3 for coordinates in spherical_coordinates) \
    and all(coordinates[2] <= 1e-3 for coordinates in spherical_coordinates):

        if len(axis_of_rotation) == principal_axis.fold:
            return True

    return False


def check_sigma_h(reflection_symmetry):
    for reflection in reflection_symmetry:
        if Vector.theta(reflection.vector) % pi <= 1e-3:
            return True
    return False


def check_icosahedron(rotation_symmetry):
    if any([vector.fold == 5 for vector in rotation_symmetry]):
        return True
    else:
        return False


def check_high_symmetry(rotation_symmetry):
    i = 0
    for rotation in rotation_symmetry:
        if rotation.fold > 2:
            i += 1
    if i > 2:
        return True
    else:
        return False


def check_inversion_symmetry(nuclei_array):

    nuclei_array_inverse = []
    for nuclei in nuclei_array:
        coordinate_inverse = (-nuclei.coordinates[0], -nuclei.coordinates[1], -nuclei.coordinates[2])
        nuclei_inverse = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinate_inverse)
        nuclei_array_inverse.append(nuclei_inverse)

    for nuclei in nuclei_array:
        for i, nuclei_inverse in enumerate(nuclei_array_inverse):
            if Vector.distance(nuclei.coordinates, nuclei_inverse.coordinates) <= 1e-3 \
            and (nuclei.charge - nuclei_inverse.charge) == 0.0:
                break
            if i == len(nuclei_array_inverse) - 1:
                return False
    return True


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


def return_principal_axis(rotation_symmetry):
    for rotation in rotation_symmetry:
        if Vector.theta(rotation.vector) % pi <= 1e-3:
            return rotation


def standard_orientation(nuclei_array, rotation_symmetry, reflection_symmetry):

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

        quaternion = quaternion_rotate_from_phi(first_highest_symmetry.vector, 0.0)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        quaternion = quaternion_rotate_to_z_axis(second_highest_symmetry.vector)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

    elif len(rotation_symmetry) == 1 and len(reflection_symmetry) >= 1:
        first_highest_symmetry = rotation_symmetry[0]
        reflection_d = None
        sigma_d = False

        quaternion = quaternion_rotate_to_z_axis(first_highest_symmetry.vector)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        for reflection in reflection_symmetry:
            if Vector.phi(reflection.vector) > 1e-3:
                reflection_d = reflection
                sigma_d = True
                break

        if sigma_d:
            quaternion = quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
            rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

    elif len(rotation_symmetry) == 0 and len(reflection_symmetry) > 1:
        reflection_h = reflection_symmetry[0]
        reflection_d = reflection_symmetry[1]

        quaternion = quaternion_rotate_to_z_axis(reflection_h.vector)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

        quaternion = quaternion_rotate_from_phi(reflection_d.vector, pi / 2)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

    elif len(rotation_symmetry) == 0 and len(reflection_symmetry) == 1:
        quaternion = quaternion_rotate_to_z_axis(reflection_symmetry[0].vector)
        rotate_all_vectors(quaternion, rotation_symmetry, reflection_symmetry, nuclei_array)

    return nuclei_array, rotation_symmetry, reflection_symmetry


def rotate_all_vectors(quaternion_i, rotation_symmetry_list_i, reflection_symmetry_list_i, nuclei_array_i):
    for rotation_i in rotation_symmetry_list_i:
        rotation_i.vector = Vector.quaternion_rotation(quaternion_i, rotation_i.vector)
    for reflection_i in reflection_symmetry_list_i:
        reflection_i.vector = Vector.quaternion_rotation(quaternion_i, reflection_i.vector)
    for nuclei_i in nuclei_array_i:
        nuclei_i.coordinates = Vector.quaternion_rotation(quaternion_i, nuclei_i.coordinates)


def quaternion_rotate_from_phi(symmetry_vector, angle):
    vector = (0.0, 0.0, 1.0)
    phi = - Vector.phi(symmetry_vector) + angle
    quaternion = Vector.create_quaternion(vector, phi)
    return quaternion


def quaternion_rotate_to_z_axis(symmetry_vector):
    vector = (- symmetry_vector[1], symmetry_vector[0], 0.0)
    theta = - Vector.theta(symmetry_vector)
    quaternion = Vector.create_quaternion(vector, theta)
    return quaternion


def brute_force_symmetry(nuclei_array):
    nuclei_array = remove_center_nuclei(nuclei_array)
    corners = vertices(nuclei_array)
    edge_center = center_two_vertices(nuclei_array)
    faces_center = center_three_vertices(nuclei_array)
    cross_corners_corners = cross_products_vertices_vertices(corners)
    cross_edge_center_corners = cross_product_center_edge_vertices(corners)

    rotation_symmetry = brute_force_rotation_symmetry(nuclei_array, corners, edge_center, faces_center,
    cross_corners_corners, cross_edge_center_corners)
    reflection_symmetry = brute_force_reflection_symmetry(nuclei_array, rotation_symmetry,
    cross_corners_corners, cross_edge_center_corners)

    return rotation_symmetry, reflection_symmetry


def remove_center_nuclei(nuclei_array):
    for i, nuclei in enumerate(nuclei_array):
        if Vector.rho(nuclei.coordinates) <= 1e-3:
            nuclei_array.pop(i)
    return nuclei_array


def vertices(nuclei_array):
    corner = []
    for nuclei in nuclei_array:
        coordinates = Vector.normalize(nuclei.coordinates)
        corner.append(coordinates)
    return corner


def center_two_vertices(nuclei_array):
    center_of_edge = []
    for nuclei_i in nuclei_array:
        for nuclei_j in nuclei_array:
            if nuclei_i is not nuclei_j:
                axis_i = nuclei_i.coordinates
                axis_j = nuclei_j.coordinates
                axis_edge = Vector.add(axis_i, axis_j)
                if Vector.rho(axis_edge) > 1e-3:
                    axis_edge = Vector.normalize(axis_edge)
                    center_of_edge.append(axis_edge)
    return center_of_edge


def center_three_vertices(nuclei_array):
    center_of_faces = []
    for nuclei_i in nuclei_array:
        for nuclei_j in nuclei_array:
            for nuclei_k in nuclei_array:
                if (nuclei_i is not nuclei_j) and (nuclei_i is not nuclei_k) and (nuclei_j is not nuclei_k):
                    axis_i = nuclei_i.coordinates
                    axis_j = nuclei_j.coordinates
                    axis_k = nuclei_k.coordinates
                    axis_face = Vector.add(Vector.add(axis_i, axis_j), axis_k)
                    if Vector.rho(axis_face) > 1e-3:
                        axis_face = Vector.normalize(axis_face)
                        center_of_faces.append(axis_face)
    return center_of_faces


def cross_products_vertices_vertices(corner):
    corner = remove_duplicate(corner)
    cross_products = []
    for axis_i in corner:
        for axis_j in corner:
            if axis_i is not axis_j:
                axis_cross = Vector.cross(axis_i, axis_j)
                if Vector.rho(axis_cross) > 1e-3:
                    axis_cross = Vector.normalize(axis_cross)
                    cross_products.append(axis_cross)
    return cross_products


def cross_product_center_edge_vertices(corner):
    corner = remove_duplicate(corner)
    cross_products = []
    for axis_i in corner:
        for axis_j in corner:
            for axis_k in corner:
                if (axis_i is not axis_j) and (axis_i is not axis_k) and (axis_j is not axis_k):
                    axis_l = Vector.add(axis_i, axis_j)
                    axis_cross = Vector.cross(axis_l, axis_k)
                    if Vector.rho(axis_cross) > 1e-3:
                        axis_cross = Vector.normalize(axis_cross)
                        cross_products.append(axis_cross)
    return cross_products


def brute_force_reflection_symmetry(nuclei_array, rotation_symmetry, cross_corners_corners,
    cross_edge_center_corners):

    # rotate all orthogonal vectors by principal axis by twice it's n-fold
    vector_cross = remove_duplicate(cross_corners_corners + cross_edge_center_corners)
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

    reflection_planes = []
    if len(vector_cross) > 0:

        total_vectors_cross = cross_corners_corners + vectors_cross_rotated
        vectors_reflection_plane = remove_duplicate(total_vectors_cross)

        # create householder matrices
        householder_matrices = []
        for planes in vectors_reflection_plane:
            planes = np.matrix(planes)
            householder_matrices.append(np.identity(3) - 2 * planes.T * planes)

        # check each reflection
        for i, matrix in enumerate(householder_matrices):
            if check_reflection(nuclei_array, matrix):
                planes_of_reflection = ReflectionSymmetry(vectors_reflection_plane[i])
                reflection_planes.append(planes_of_reflection)

    return reflection_planes


def check_reflection(nuclei_array_i, householder_matrix):
    nuclei_array_copy = []
    for nuclei in nuclei_array_i:
        coordinates = householder_matrix * np.matrix(nuclei.coordinates).T
        coordinates = tuple(coordinates.T.tolist()[0])
        nuclei_copy = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinates)
        nuclei_array_copy.append(nuclei_copy)
    for nuclei_i in nuclei_array_i:
        for j, nuclei_j in enumerate(nuclei_array_copy):
            if Vector.distance(nuclei_i.coordinates, nuclei_j.coordinates) <= 1e-3 \
            and (nuclei_i.charge - nuclei_j.charge) == 0.0:
                break
            if j == len(nuclei_array_copy) - 1:
                return False
    return True


def brute_force_rotation_symmetry(nuclei_array, corner, edge_center, faces_center, cross_corners_corners,
    cross_edge_corners_vertices):

    axis_of_rotations_i = remove_duplicate(corner + edge_center + faces_center + cross_corners_corners
    + cross_edge_corners_vertices)
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
                if check_quaternion(nuclei_array, quaternion_matrix[i, j]):
                    n_fold_symmetry_i[j] = i + 2

        # create the rotation symmetry object and return them if the symmetry if > 1-fold
        for i, symmetry in enumerate(n_fold_symmetry_i):
            if symmetry != 1:
                rotation_symmetry = RotationSymmetry(n_fold_symmetry_i[i], axis_of_rotations_i[i])
                axis_of_rotations_j.append(rotation_symmetry)

    return axis_of_rotations_j


def check_quaternion(nuclei_array_i, quaternion):
    nuclei_array_copy = []
    for nuclei in nuclei_array_i:
        coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
        nuclei_copy = Nuclei(nuclei.element, nuclei.charge, nuclei.mass, coordinates)
        nuclei_array_copy.append(nuclei_copy)
    for nuclei_i in nuclei_array_i:
        for k, nuclei_k in enumerate(nuclei_array_copy):
            if Vector.distance(nuclei_i.coordinates, nuclei_k.coordinates) <= 1e-3 \
            and (nuclei_i.charge - nuclei_k.charge) == 0.0:
                break
            if k == len(nuclei_array_copy) - 1:
                return False
    return True


def remove_duplicate(axis_of_rotations):
    axis_of_rotations_j = []
    for axis_j in axis_of_rotations:
        if check_array(axis_j, axis_of_rotations_j):
            axis_of_rotations_j.append(axis_j)
    return axis_of_rotations_j


def check_array(axis, axis_of_rotations_i):
    for axis_i in axis_of_rotations_i:
        if Vector.rho(Vector.cross(axis, axis_i)) <= 1e-3:
            return False
    return True


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
