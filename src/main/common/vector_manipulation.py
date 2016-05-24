from math import sqrt
from math import cos, sin, acos, atan2


class Vector:

    @staticmethod
    def add(r_1, r_2):
        return r_1[0] + r_2[0], r_1[1] + r_2[1], r_1[2] + r_2[2]

    @staticmethod
    def minus(r_1, r_2):
        return r_1[0] - r_2[0], r_1[1] - r_2[1], r_1[2] - r_2[2]

    @staticmethod
    def distance(r_1, r_2):
        r_12 = sqrt((r_1[0] - r_2[0])**2 + (r_1[1] - r_2[1])**2 + (r_1[2] - r_2[2])**2)
        return r_12

    @staticmethod
    def cross(r_1, r_2):
        i = r_1[1] * r_2[2] - r_1[2] * r_2[1]
        j = r_1[2] * r_2[0] - r_1[0] * r_2[2]
        k = r_1[0] * r_2[1] - r_1[1] * r_2[0]
        return i, j, k

    @staticmethod
    def quaternion_multi(q_1, q_2):
        a_1, b_1, c_1, d_1 = q_1
        a_2, b_2, c_2, d_2 = q_2

        a = (a_1 * a_2) - (b_1 * b_2) - (c_1 * c_2) - (d_1 * d_2)
        b = (a_1 * b_2) + (b_1 * a_2) + (c_1 * d_2) - (d_1 * c_2)
        c = (a_1 * c_2) - (b_1 * d_2) + (c_1 * a_2) + (d_1 * b_2)
        d = (a_1 * d_2) + (b_1 * c_2) - (c_1 * b_2) + (d_1 * a_2)

        return a, b, c, d

    @staticmethod
    def quaternion_conj(q):
        a, b, c, d = q
        return a, -b, -c, -d

    @classmethod
    def create_quaternion(cls, axis, angle):
        if all(abs(value) <= 1e-3 for value in axis):
            return 1.0, 0.0, 0.0, 0.0
        else:
            axis = cls.normalize(axis)
            a = cos(angle / 2)
            b = axis[0] * sin(angle / 2)
            c = axis[1] * sin(angle / 2)
            d = axis[2] * sin(angle / 2)
        return a, b, c, d

    @classmethod
    def quaternion_rotation(cls, quaternion, point):
        point = (0.0, point[0], point[1], point[2])
        point = cls.quaternion_multi(cls.quaternion_multi(quaternion, point), cls.quaternion_conj(quaternion))
        return point[1], point[2], point[3]

    @staticmethod
    def normalize(r):
        magnitude = sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        r_x = r[0] / magnitude
        r_y = r[1] / magnitude
        r_z = r[2] / magnitude
        return r_x, r_y, r_z

    @staticmethod
    def rho(r):
        return sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)

    @staticmethod
    def theta(r):
        rho = sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
        if rho <= 1e-3:
            return 0.0
        else:
            return acos(r[2] / rho)

    @staticmethod
    def phi(r):
        if abs(r[1]) <= 1e-3:
            return 0.0
        else:
            return atan2(r[1], r[0])

    @classmethod
    def cartesian_to_spherical(cls, r):
        rho = cls.rho(r)
        theta = cls.theta(r)
        phi = cls.phi(r)
        return rho, theta, phi

    @staticmethod
    def gaussian(a, r_1, b, r_2):
        i = (a * r_1[0] + b * r_2[0]) / (a + b)
        j = (a * r_1[1] + b * r_2[1]) / (a + b)
        k = (a * r_1[2] + b * r_2[2]) / (a + b)
        return i, j, k