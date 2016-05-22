from math import isclose
from math import pi
from src.main.common.vector_manipulation import Vector


class Symmetry:

    @classmethod
    def check_sym(cls, basis1, basis2, basis3, basis4):
        r_1 = basis1.coordinates
        r_2 = basis2.coordinates
        r_3 = basis3.coordinates
        r_4 = basis4.coordinates
        l_1 = basis1.integral_exponents
        l_2 = basis2.integral_exponents
        l_3 = basis3.integral_exponents
        l_4 = basis4.integral_exponents

        if cls.even_odd(l_1[0]) * cls.even_odd(l_2[0]) * cls.even_odd(l_3[0]) * cls.even_odd(l_4[0]) == -1:
            if isclose(r_1[0], r_2[0], abs_tol=1e-3) and isclose(r_1[0], r_3[0], abs_tol=1e-3) \
            and isclose(r_1[0], r_4[0], abs_tol=1e-3):
                return False
            else:
                return True
        elif cls.even_odd(l_1[1]) * cls.even_odd(l_2[1]) * cls.even_odd(l_3[1]) * cls.even_odd(l_4[1]) == -1:
            if isclose(r_1[1], r_2[1], abs_tol=1e-3) and isclose(r_1[1], r_3[1], abs_tol=1e-3) \
            and isclose(r_1[1], r_4[1], abs_tol=1e-3):
                return False
            else:
                return True
        elif cls.even_odd(l_1[2]) * cls.even_odd(l_2[2]) * cls.even_odd(l_3[2]) * cls.even_odd(l_4[2]) == -1:
            if isclose(r_1[2], r_2[2], abs_tol=1e-3) and isclose(r_1[2], r_3[2], abs_tol=1e-3) \
            and isclose(r_1[2], r_4[2], abs_tol=1e-3):
                return False
            else:
                return True
        else:
            return True

    @staticmethod
    def even_odd(num):
        if num % 2 == 0:
            return 1
        else:
            return -1

    @staticmethod
    def sort(i, j, k, l):
        a = i
        b = j
        c = k
        d = l
        if a > b:
            a, b = b, a
        if c > d:
            c, d = d, c
        if a > c:
            a, c = c, a
            b, d = d, b
        if a == c:
            if b > d:
                a, c = c, a
                b, d = d, b
        return a, b, c, d

    @classmethod
    def point_group(cls, nuclei_array):
        number_of_nuclei = len(nuclei_array)

        # if only one nucleus then center the nuclei and skip symmetry checks
        if number_of_nuclei == 1:
            nuclei = nuclei_array[0]
            nuclei.coordinates = (0.0, 0.0, 0.0)
            return [nuclei]

        else:
            # put the molecule in the standard orientation
            nuclei_array = cls.standard_orientation(nuclei_array)

            # determine the point group first check if molecule is linear
            if cls.check_linear(nuclei_array):

                # check for inversion symmetry
                if cls.check_inversion_symmetry(nuclei_array):
                    print('true')
                else:
                    print('false')

            # check if there are two or more c_n n > 2
            else:
                pass

        return nuclei_array

    @staticmethod
    def standard_orientation(nuclei_array):
        number_of_nuclei = len(nuclei_array)
        nuclei_array_1 = []
        nuclei_array_2 = []
        center = [0, 0, 0]
        quaternion = (1.0, 0.0, 0.0, 0.0)

        # find the center of the molecule
        for nuclei in nuclei_array:
            center[0] += nuclei.coordinates[0] / number_of_nuclei
            center[1] += nuclei.coordinates[1] / number_of_nuclei
            center[2] += nuclei.coordinates[2] / number_of_nuclei

        # move the molecule so that the center is at the origin
        for nuclei in nuclei_array:
            x = nuclei.coordinates[0] - center[0]
            y = nuclei.coordinates[1] - center[1]
            z = nuclei.coordinates[2] - center[2]
            nuclei.coordinates = (x, y, z)
            nuclei_array_1.append(nuclei)

        # sort array by distance from origin
        nuclei_array_1.sort(key=lambda nucleus: Vector.rho(nucleus.coordinates))

        # check if a nuclei is close to the origin
        if Vector.rho(nuclei_array_1[0].coordinates) <= 1e-3:

            # move the closest nuclei to the origin
            for i, nuclei in enumerate(nuclei_array_1):
                if i == 0:
                    translation = nuclei.coordinates
                    nuclei.coordinates = (0.0, 0.0, 0.0)
                    nuclei_array_2.append(nuclei)
                else:
                    x = nuclei.coordinates[0] - translation[0]
                    y = nuclei.coordinates[1] - translation[1]
                    z = nuclei.coordinates[2] - translation[2]
                    nuclei.coordinates = (x, y, z)
                    nuclei_array_2.append(nuclei)

            # sort by the angle theta keeping the nuclei at the origin at index 0
            initial_nuclei = nuclei_array_2[0]
            nuclei_array_2 = nuclei_array_2[1:]
            nuclei_array_2.sort(key=lambda nucleus: Vector.theta(nucleus.coordinates))
            nuclei_array_2.insert(0, initial_nuclei)
            nuclei_array_1 = []

            # rotate the molecule so that the nuclei at index 1 is on the z-axis
            for i, nuclei in enumerate(nuclei_array_2):
                if i == 0:
                    nuclei_array_1.append(nuclei)
                elif i == 1:
                    vector = (- nuclei.coordinates[1], nuclei.coordinates[0], 0.0)
                    theta = - Vector.theta(nuclei.coordinates)
                    if all(abs(value) <= 1e-3 for value in vector):
                        quaternion = (1.0, 0.0, 0.0, 0.0)
                    else:
                        quaternion = Vector.create_quaternion(vector, theta)
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_1.append(nuclei)
                else:
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_1.append(nuclei)

            # sort by the angle phi keeping nuclei at origin and nuclei on z-axis at index 0 and 1
            initial_nuclei = [nuclei_array_1[0], nuclei_array_1[1]]
            nuclei_array_1 = nuclei_array_1[2:]
            nuclei_array_1.sort(key=lambda nucleus: abs(Vector.phi(nucleus.coordinates)))
            nuclei_array_1 = initial_nuclei + nuclei_array_1
            nuclei_array_2 = []

            # rotate the molecule so that the nucleus at index 2 is on the xz-plane
            for i, nuclei in enumerate(nuclei_array_1):
                if i == 0 or i == 1:
                    nuclei_array_2.append(nuclei)
                elif i == 2:
                    vector = (0.0, 0.0, 1.0)
                    phi = - Vector.phi(nuclei.coordinates)
                    quaternion = Vector.create_quaternion(vector, phi)
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_2.append(nuclei)
                else:
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_2.append(nuclei)

            return nuclei_array_2

        # if there is not a nuclei close to the origin
        elif not Vector.rho(nuclei_array_1[0].coordinates) <= 1e-3:
            nuclei_array_2 = nuclei_array_1
            nuclei_array_1 = []

            # rotate the molecule so that the nuclei at index 0 is on the z-axis
            for i, nuclei in enumerate(nuclei_array_2):
                if i == 0:
                    vector = (- nuclei.coordinates[1], nuclei.coordinates[0], 0.0)
                    theta = - Vector.theta(nuclei.coordinates)
                    if all(abs(value) <= 1e-3 for value in vector):
                        quaternion = (1.0, 0.0, 0.0, 0.0)
                    else:
                        quaternion = Vector.create_quaternion(vector, theta)
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_1.append(nuclei)
                else:
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_1.append(nuclei)

            # sort by angle phi keeping the nucleus on the z-axis at index 0
            initial_nuclei = nuclei_array_1[0]
            nuclei_array_1 = nuclei_array_1[1:]
            nuclei_array_1.sort(key=lambda nucleus: abs(Vector.phi(nucleus.coordinates)))
            nuclei_array_1.insert(0, initial_nuclei)
            nuclei_array_2 = []

            # rotate the molecule so that the nucleus at index 1 is on the xz-plane
            for i, nuclei in enumerate(nuclei_array_1):
                if i == 0:
                    nuclei_array_2.append(nuclei)
                elif i == 1:
                    vector = (0.0, 0.0, 1.0)
                    phi = - Vector.phi(nuclei.coordinates)
                    quaternion = Vector.create_quaternion(vector, phi)
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_2.append(nuclei)
                else:
                    nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                    nuclei_array_2.append(nuclei)

            return nuclei_array_2

    @staticmethod
    def check_linear(nuclei_array):
        spherical_coordinates = []

        for nuclei in nuclei_array:
            nuclei.coordinates = Vector.cartesian_to_spherical(nuclei.coordinates)
            spherical_coordinates.append(nuclei)

        if all(nuclei.coordinates[1] % pi <= 1e-3 for nuclei in spherical_coordinates) \
        and all(nuclei.coordinates[2] <= 1e-3 for nuclei in spherical_coordinates):
            return True
        else:
            return False

    @staticmethod
    def check_inversion_symmetry(nuclei_array):

        for i, nuclei in enumerate(nuclei_array):
            coordinate_inverse = (- nuclei.coordinates[0], - nuclei.coordinates[1], - nuclei.coordinates[2])
            check_list = nuclei_array
            check_list.pop(i)

            if not any(nucleus.element == nuclei.element for nucleus in check_list) \
            or not any(nucleus.coordinates == coordinate_inverse for nucleus in check_list):
                return False

        return True
