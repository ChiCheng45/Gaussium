from src.main.common import rho, theta, phi
from src.main.common import cartesian_to_spherical
from src.main.common import create_quaternion
from src.main.common import quaternion_rotation
from src.main.common import quaternion_multi
from src.main.factory import SymmetryFactory
from src.main.objects import PointGroup
from src.main.objects import Oh, D4h, C4v
from math import pi
import heapq


class MoleculeFactory:

    def __init__(self, symmetry=False, error=1e-2):
        self.symmetry = symmetry
        self.error = error
        self.symmetry_factory = SymmetryFactory(error)

    def create(self, nuclei_array):

        if not self.symmetry:
            return nuclei_array, PointGroup([], [], [], [], 'C_{1}')

        nuclei_array = self.center_molecule(nuclei_array)

        if len(nuclei_array) == 1:                                      # Point
            return nuclei_array, Oh()

        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = self.standard_orientation(nuclei_array, rotation, reflection)

        if self.check_linear(nuclei_array):                             # Linear
            if len(inversion) == 1:
                return nuclei_array, D4h()
            else:
                return nuclei_array, C4v()

        if self.check_high_symmetry(rotation):                          # Polyhedral
            if not len(inversion) == 1:
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'T_{d}')
            elif any([vector.fold == 5 for vector in rotation]):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'I_{h}')
            else:
                return nuclei_array, PointGroup(rotation, reflection, inversion, improper, 'O_{h}')

        if len(rotation) == 0:                                          # Nonaxial
            if self.check_sigma_h(reflection):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{s}')
            elif len(inversion) == 1:
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{i}')
            else:
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{1}')

        n = self.return_principal_axis(rotation).fold

        if self.check_n_two_fold_perpendicular_to_n_fold(rotation):     # Dihedral
            if self.check_sigma_h(reflection):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'D_{'+str(n)+'h}')
            elif self.check_n_sigma_v(n, reflection):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'D_{'+str(n)+'d}')
            else:
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'D_{'+str(n)+'}')

        else:                                                           # Cyclic
            if self.check_sigma_h(reflection):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{'+str(n)+'h}')
            elif self.check_n_sigma_v(n, reflection):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{'+str(n)+'v}')
            elif self.check_improper_rotation(n, improper):
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'S_{2'+str(n)+'}')
            else:
                return nuclei_array, PointGroup(rotation, reflection, improper, inversion, 'C_{'+str(n)+'}')

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

    def standard_orientation(self, nuclei_array, rotation_symmetry, reflection_symmetry):
        vector_i = vector_j = (0.0, 0.0, 0.0)

        if len(rotation_symmetry) > 1:
            highest_n_folds = heapq.nlargest(2, [rotation.fold for rotation in rotation_symmetry])
            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[0]:
                    vector_i = rotation.vector
                    break
            for rotation in rotation_symmetry:
                if rotation.fold == highest_n_folds[1] and rotation.vector != vector_i:
                    vector_j = rotation.vector
                    break

        if len(rotation_symmetry) == 1:
            vector_i = rotation_symmetry[0].vector
            for reflection in reflection_symmetry:
                if phi(reflection.vector) > self.error:
                    vector_j = reflection.vector
                    break

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) > 1:
            vector_i = reflection_symmetry[0].vector
            vector_j = reflection_symmetry[1].vector

        if rho(nuclei_array[0].coordinates) > 1e-3:
            i = 0
        else:
            i = 1

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) == 1:
            vector_i = reflection_symmetry[0].vector
            vector_j = nuclei_array[i].coordinates

        if len(rotation_symmetry) == 0 and len(reflection_symmetry) == 0:
            vector_i = nuclei_array[i].coordinates
            vector_j = nuclei_array[i + 1].coordinates

        if rho(vector_i) <= self.error:
            vector_i = (1.0, 0.0, 0.0)

        quaternion_i = create_quaternion((-vector_i[1], vector_i[0], 0.0), -theta(vector_i))
        quaternion_j = create_quaternion((0.0, 0.0, 1.0), -phi(vector_j))
        quaternion = quaternion_multi(quaternion_j, quaternion_i)

        for rotation in rotation_symmetry:
            rotation.vector = quaternion_rotation(quaternion, rotation.vector)
        for reflection in reflection_symmetry:
            reflection.vector = quaternion_rotation(quaternion, reflection.vector)
        for nuclei in nuclei_array:
            nuclei.coordinates = quaternion_rotation(quaternion, nuclei.coordinates)

        return nuclei_array, rotation_symmetry, reflection_symmetry

    def check_improper_rotation(self, n, improper_rotations):
        for improper_rotation in improper_rotations:
            if improper_rotation.fold == 2 * n:
                return True
        return False

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
