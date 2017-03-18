import heapq, os
from contextlib import redirect_stdout
import numpy as np
from src.common import create_quaternion
from src.common import normalize
from src.common import quaternion_rotation
from src.common import read_basis_set_file
from src.common import theta, phi


class NelderMead:

    def __init__(self, basis_file, energy_object, nuclei_list, tau=0.25, threshold=1e-6):
        self.basis_file = basis_file
        self.energy_object = energy_object
        self.first_nuclei, self.nuclei_list = self.reorient_molecule(nuclei_list)
        if len(self.nuclei_list) == 1:
            self.m = 1
        else:
            self.m = len(self.nuclei_list) * 3 - 3
        self.alpha = 1
        self.beta = 2 + (2 / self.m)
        self.gamma = 0.75 + (1 / (2 * self.m))
        self.delta = 0.5 - (1 / self.m)
        self.tau = tau
        self.threshold = threshold

    def optimize(self):
        print('\n*************************************************************************************************')
        print('\nELECTRONIC STRUCTURE METHOD: {}'.format(self.energy_object.method))
        print('\nBEGIN NELDER-MEAD OPTIMIZATION')

        simplex_matrix, energy_list = self.build_initial_simplex()

        while True:
            max_energy_m, max_energy_l = heapq.nlargest(2, energy_list)
            min_energy = min(energy_list)
            max_energy_index_m = energy_list.index(max_energy_m)

            max_energy_points = np.matrix(simplex_matrix[:, max_energy_index_m])
            centroid = np.matrix((1 / self.m) * (np.sum(simplex_matrix, axis=1) - max_energy_points))
            reflection = centroid + self.alpha * (centroid - max_energy_points)
            reflection_energy = self.calculate_energy(reflection)
            if min_energy <= reflection_energy < max_energy_l:
                simplex_matrix, energy_list = self.replace_max_energy(
                    simplex_matrix, energy_list, max_energy_index_m, reflection, reflection_energy
                )

            if reflection_energy < min_energy:
                expansion = centroid + self.beta * (reflection - centroid)
                expansion_energy = self.calculate_energy(expansion)

                if expansion_energy < reflection_energy:
                    simplex_matrix, energy_list = self.replace_max_energy(
                        simplex_matrix, energy_list, max_energy_index_m, expansion, expansion_energy
                    )
                else:
                    simplex_matrix, energy_list = self.replace_max_energy(
                        simplex_matrix, energy_list, max_energy_index_m, reflection, reflection_energy
                    )

            if max_energy_l <= reflection_energy < max_energy_m:
                outside_contraction = centroid + self.gamma * (reflection - centroid)
                outside_energy = self.calculate_energy(outside_contraction)

                if outside_energy <= reflection_energy:
                    simplex_matrix, energy_list = self.replace_max_energy(
                        simplex_matrix, energy_list, max_energy_index_m, outside_contraction, outside_energy
                    )
                else:
                    simplex_matrix = self.shrink_simplex(simplex_matrix, energy_list)
                    energy_list = self.calculate_simplex_energies(simplex_matrix)

            if max_energy_m <= reflection_energy:
                inside_contraction = centroid - self.gamma * (reflection - centroid)
                inside_energy = self.calculate_energy(inside_contraction)

                if inside_energy < max_energy_m:
                    simplex_matrix, energy_list = self.replace_max_energy(
                        simplex_matrix, energy_list, max_energy_index_m, inside_contraction, inside_energy
                    )
                else:
                    simplex_matrix = self.shrink_simplex(simplex_matrix, energy_list)
                    energy_list = self.calculate_simplex_energies(simplex_matrix)

            print("ENERGY: {}, DEV: {}".format(min(energy_list), np.std(energy_list)))

            if np.std(energy_list) <= self.threshold:
                break

        return min(energy_list)

    def build_initial_simplex(self):

        initial_points = []
        for i, nuclei in enumerate(self.nuclei_list):
            if i == 0:
                initial_points.append(nuclei.coordinates[2])
            elif i == 1:
                initial_points.append(nuclei.coordinates[1])
                initial_points.append(nuclei.coordinates[2])
            else:
                initial_points.append(nuclei.coordinates[0])
                initial_points.append(nuclei.coordinates[1])
                initial_points.append(nuclei.coordinates[2])
        initial_points = np.matrix(initial_points).T

        simplex_matrix = np.copy(initial_points)
        for i in range(self.m):
            points = np.copy(initial_points)
            value = points.item(i, 0)
            points.itemset((i, 0), value + self.tau)
            simplex_matrix = np.concatenate((simplex_matrix, points), axis=1)
        simplex_matrix = np.matrix(simplex_matrix)

        energy_list = self.calculate_simplex_energies(simplex_matrix)
        return simplex_matrix, energy_list

    def calculate_simplex_energies(self, simplex_matrix):
        energy_list = []
        for i in range(self.m + 1):
            energy = self.calculate_energy(simplex_matrix[:, i])
            energy_list.append(energy)
        return energy_list

    def calculate_energy(self, coordinate_list):

        j = 0
        for i, nuclei in enumerate(self.nuclei_list):
            if i == 0:
                nuclei.coordinates = (0.0, 0.0, coordinate_list.item(j))
                j += 1
            elif i == 1:
                nuclei.coordinates = (0.0, coordinate_list.item(j), coordinate_list.item(j+1))
                j += 2
            else:
                nuclei.coordinates = (coordinate_list.item(j), coordinate_list.item(j+1), coordinate_list.item(j+2))
                j += 3

        with redirect_stdout(open(os.devnull, "w")):
            basis_set = read_basis_set_file(self.basis_file, self.first_nuclei + self.nuclei_list)
            energy = self.energy_object.calculate_energy(self.first_nuclei + self.nuclei_list, basis_set)

        return energy

    def replace_max_energy(self, simplex_matrix, energy_list, max_energy_index, new_points, new_energy):
        del energy_list[max_energy_index]
        simplex_matrix = np.delete(simplex_matrix, max_energy_index, axis=1)
        energy_list.append(new_energy)
        simplex_matrix = np.concatenate((simplex_matrix, new_points), axis=1)
        return simplex_matrix, energy_list

    def shrink_simplex(self, simplex_matrix, energy_list):
        min_energy_index = energy_list.index(min(energy_list))
        min_points = np.matrix(simplex_matrix[:, min_energy_index])
        min_points_matrix = np.tile(min_points, (self.m + 1, ))
        simplex_matrix = min_points_matrix + self.delta * (simplex_matrix - min_points_matrix)
        return simplex_matrix

    def reorient_molecule(self, nuclei_list):
        first_nuclei = nuclei_list.pop(0)
        coordinates = first_nuclei.coordinates
        first_nuclei.coordinates = (0, 0, 0)

        for nuclei in nuclei_list:
            x = nuclei.coordinates[0] - coordinates[0]
            y = nuclei.coordinates[1] - coordinates[1]
            z = nuclei.coordinates[2] - coordinates[2]
            nuclei.coordinates = (x, y, z)

        if len(nuclei_list) >= 1:
            second_nuclei = nuclei_list[0]
            coordinates = normalize(second_nuclei.coordinates)

            quaternion = create_quaternion((-coordinates[1], coordinates[0], 0.0), -theta(coordinates))
            for nuclei in nuclei_list:
                nuclei.coordinates = quaternion_rotation(quaternion, nuclei.coordinates)

        if len(nuclei_list) >= 2:
            third_nuclei = nuclei_list[1]
            coordinates = normalize(third_nuclei.coordinates)

            quaternion = create_quaternion((0.0, 0.0, 1.0), -phi(coordinates) + np.pi/2)
            for nuclei in nuclei_list:
                nuclei.coordinates = quaternion_rotation(quaternion, nuclei.coordinates)

        return [first_nuclei], nuclei_list
