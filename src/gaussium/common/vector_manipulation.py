from math import sqrt, cos, sin, acos, atan2
import numpy as np


def vector_add(r_1, r_2):
    return r_1[0] + r_2[0], r_1[1] + r_2[1], r_1[2] + r_2[2]


def vector_minus(r_1, r_2):
    return r_1[0] - r_2[0], r_1[1] - r_2[1], r_1[2] - r_2[2]


def vector_multi(r, x):
    return r[0] * x, r[1] * x, r[2] * x


def vector_divide(r, x):
    return r[0] / x, r[1] / x, r[2] / x


def coordinate_distance(r_1, r_2):
    return sqrt((r_1[0] - r_2[0])**2 + (r_1[1] - r_2[1])**2 + (r_1[2] - r_2[2])**2)


def dot_product(r_1, r_2):
    i = r_1[0] * r_2[0]
    j = r_1[1] * r_2[1]
    k = r_1[2] * r_2[2]
    return i + j + k


def cross_product(r_1, r_2):
    i = r_1[1] * r_2[2] - r_1[2] * r_2[1]
    j = r_1[2] * r_2[0] - r_1[0] * r_2[2]
    k = r_1[0] * r_2[1] - r_1[1] * r_2[0]
    return i, j, k


def create_quaternion(axis, angle):
    if all(abs(value) <= 1e-3 for value in axis):
        return 1.0, 0.0, 0.0, 0.0
    else:
        axis = normalize(axis)
        a = cos(angle / 2)
        b = axis[0] * sin(angle / 2)
        c = axis[1] * sin(angle / 2)
        d = axis[2] * sin(angle / 2)
    return a, b, c, d


def quaternion_rotation(quaternion, point):
    point = (0.0, point[0], point[1], point[2])
    point = quaternion_multi(quaternion_multi(quaternion, point), quaternion_conjugate(quaternion))
    return point[1], point[2], point[3]


def quaternion_multi(q_1, q_2):
    a_1, b_1, c_1, d_1 = q_1
    a_2, b_2, c_2, d_2 = q_2

    a = (a_1 * a_2) - (b_1 * b_2) - (c_1 * c_2) - (d_1 * d_2)
    b = (a_1 * b_2) + (b_1 * a_2) + (c_1 * d_2) - (d_1 * c_2)
    c = (a_1 * c_2) - (b_1 * d_2) + (c_1 * a_2) + (d_1 * b_2)
    d = (a_1 * d_2) + (b_1 * c_2) - (c_1 * b_2) + (d_1 * a_2)

    return a, b, c, d


def create_householder_matrix(vector):
    planes = np.array([vector])
    return np.identity(3) - 2 * (planes.T @ planes)


def householder_matrix_reflection(coordinates, householder_matrix):
    coordinates = householder_matrix @ np.array(coordinates).T
    coordinates = tuple(coordinates.T.tolist())
    return coordinates


def quaternion_conjugate(q):
    a, b, c, d = q
    return a, -b, -c, -d


def normalize(r):
    magnitude = rho(r)
    r_x = r[0] / magnitude
    r_y = r[1] / magnitude
    r_z = r[2] / magnitude
    return r_x, r_y, r_z


def cartesian_to_spherical(r):
    return rho(r), theta(r), phi(r)


def rho(r):
    return sqrt(r[0]**2 + r[1]**2 + r[2]**2)


def theta(r):
    rho_i = rho(r)
    if rho_i <= 1e-3:
        return 0.0
    else:
        return acos(r[2] / rho_i)


def phi(r):
    if abs(r[1]) <= 1e-3:
        return 0.0
    else:
        return atan2(r[1], r[0])


def gaussian_product_coordinate(a, r_1, b, r_2):
    i = (a * r_1[0] + b * r_2[0]) / (a + b)
    j = (a * r_1[1] + b * r_2[1]) / (a + b)
    k = (a * r_1[2] + b * r_2[2]) / (a + b)
    return i, j, k
