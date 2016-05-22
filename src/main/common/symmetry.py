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

    @staticmethod
    def standard_orientation(nuclei_array):
        standard_nuclei_array = []
        translated_nuclei_1 = []
        translated_nuclei_2 = []
        number_of_nuclei = len(nuclei_array)
        quaternion = (1.0, 0.0, 0.0, 0.0)

        if number_of_nuclei == 1:
            nuclei = nuclei_array[0]
            nuclei.coordinates = (0.0, 0.0, 0.0)
            standard_nuclei_array.append(nuclei)
            return standard_nuclei_array

        else:
            # find the center of the molecule
            center = [0, 0, 0]
            for nuclei in nuclei_array:
                center[0] += nuclei.coordinates[0] / number_of_nuclei
                center[1] += nuclei.coordinates[1] / number_of_nuclei
                center[2] += nuclei.coordinates[2] / number_of_nuclei

            # move the molecule so that the center is at the origin and place into temp array 1
            for nuclei in nuclei_array:
                x = nuclei.coordinates[0] - center[0]
                y = nuclei.coordinates[1] - center[1]
                z = nuclei.coordinates[2] - center[2]
                nuclei.coordinates = (x, y, z)
                translated_nuclei_1.append(nuclei)

            # sort array by distance from origin
            translated_nuclei_1.sort(key=lambda nucleus: Vector.rho(nucleus.coordinates))

            # check if a nuclei is close to the origin
            if Vector.rho(translated_nuclei_1[0].coordinates) <= 1e-3:

                # move the closest nuclei to the origin if it is within 1e-3
                for i, nuclei in enumerate(translated_nuclei_1):
                    if i == 0:
                        translation = nuclei.coordinates
                        nuclei.coordinates = (0.0, 0.0, 0.0)
                        translated_nuclei_2.append(nuclei)
                    else:
                        x = nuclei.coordinates[0] - translation[0]
                        y = nuclei.coordinates[1] - translation[1]
                        z = nuclei.coordinates[2] - translation[2]
                        nuclei.coordinates = (x, y, z)
                        translated_nuclei_2.append(nuclei)

                # sort by the angle theta keeping the nuclei at the origin at index 0
                initial_nuclei = translated_nuclei_2[0]
                translated_nuclei_2 = translated_nuclei_2[1:]
                translated_nuclei_2.sort(key=lambda nucleus: Vector.theta(nucleus.coordinates))
                translated_nuclei_2.insert(0, initial_nuclei)
                translated_nuclei_1 = []

                # rotate the molecule so that the nuclei at index 1 is on the z-axis
                for i, nuclei in enumerate(translated_nuclei_2):
                    if i == 0:
                        translated_nuclei_1.append(nuclei)
                    elif i == 1:
                        vector = (- nuclei.coordinates[1], nuclei.coordinates[0], 0.0)
                        theta = - Vector.theta(nuclei.coordinates)
                        if all(value <= 1e-3 for value in vector):
                            quaternion = (1.0, 0.0, 0.0, 0.0)
                        else:
                            quaternion = Vector.create_quaternion(vector, theta)
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_1.append(nuclei)
                    else:
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_1.append(nuclei)

                # sort by the angle phi keeping nuclei at orgin and nuclei on z-axis at index 0 and 1
                initial_nuclei = [translated_nuclei_1[0], translated_nuclei_1[1]]
                translated_nuclei_1 = translated_nuclei_1[2:]
                translated_nuclei_1.sort(key=lambda nucleus: abs(Vector.phi(nucleus.coordinates)))
                translated_nuclei_1 = initial_nuclei + translated_nuclei_1
                translated_nuclei_2 = []

                # rotate the molecule so that the nucleus at index 2 is on the xz-plane
                for i, nuclei in enumerate(translated_nuclei_1):
                    if i == 0 or i == 1:
                        translated_nuclei_2.append(nuclei)
                    elif i == 2:
                        vector = (0.0, 0.0, 1.0)
                        phi = - Vector.phi(nuclei.coordinates)
                        quaternion = Vector.create_quaternion(vector, phi)
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_2.append(nuclei)
                    else:
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_2.append(nuclei)

            elif not Vector.rho(translated_nuclei_1[0].coordinates) <= 1e-3:
                translated_nuclei_2 = translated_nuclei_1
                translated_nuclei_1 = []

                # rotate the molecule so that the nuclei at index 0 is on the z-axis
                for i, nuclei in enumerate(translated_nuclei_2):
                    if i == 0:
                        vector = (- nuclei.coordinates[1], nuclei.coordinates[0], 0.0)
                        theta = - Vector.theta(nuclei.coordinates)
                        if all(value <= 1e-3 for value in vector):
                            quaternion = (1.0, 0.0, 0.0, 0.0)
                        else:
                            quaternion = Vector.create_quaternion(vector, theta)
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_1.append(nuclei)
                    else:
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_1.append(nuclei)

                # sort by angle phi keeping the nucleus on the z-axis at index 0
                initial_nuclei = translated_nuclei_1[0]
                translated_nuclei_1 = translated_nuclei_1[1:]
                translated_nuclei_1.sort(key=lambda nucleus: abs(Vector.phi(nucleus.coordinates)))
                translated_nuclei_1.insert(0, initial_nuclei)
                translated_nuclei_2 = []

                # rotate the molecule so that the nucleus at index 1 is on the xz-plane
                for i, nuclei in enumerate(translated_nuclei_1):
                    if i == 0:
                        translated_nuclei_2.append(nuclei)
                    elif i == 1:
                        vector = (0.0, 0.0, 1.0)
                        phi = - Vector.phi(nuclei.coordinates)
                        quaternion = Vector.create_quaternion(vector, phi)
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_2.append(nuclei)
                    else:
                        nuclei.coordinates = Vector.quaternion_rotation(quaternion, nuclei.coordinates)
                        translated_nuclei_2.append(nuclei)

            # determine the point group first change to spherical coordinates
            translated_nuclei_1 = []
            for nuclei in translated_nuclei_2:
                nuclei.coordinates = Vector.cartesian_to_spherical(nuclei.coordinates)
                translated_nuclei_1.append(nuclei)
                print(nuclei.coordinates)

            # check if molecule is linear
            if all(nuclei.coordinates[1] % pi <= 1e-3 for nuclei in translated_nuclei_1) \
            and all(nuclei.coordinates[2] <= 1e-3 for nuclei in translated_nuclei_1):

                # check if the molecule has a center of inversion
                pass

            # two or more c_n n > 2
            else:
                pass

        return standard_nuclei_array
