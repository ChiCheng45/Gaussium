from src.main.common import read_basis_set_file
from contextlib import redirect_stdout
import numpy as np
import itertools, heapq, os


class NelderMead:

    def __init__(self, basis_file, energy_object, nuclei_list, tau=0.25, threshold=1e-3):
        self.basis_file = basis_file
        self.energy_object = energy_object
        self.nuclei_list = nuclei_list
        self.m = len(nuclei_list) * 3
        self.alpha = 1
        self.beta = 2 + (2 / self.m)
        self.gamma = 0.75 + (1 / (2 * self.m))
        self.delta = 0.5 - (1 / self.m)
        self.tau = tau
        self.threshold = threshold

    def optimize(self):
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

            print(min(energy_list), np.std(energy_list))

            if np.std(energy_list) <= self.threshold:
                break

        return min(energy_list)

    def build_initial_simplex(self):
        initial_points = np.matrix(
            list(itertools.chain.from_iterable([nuclei.coordinates for nuclei in self.nuclei_list]))
        ).T

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
        for i, nuclei in enumerate(self.nuclei_list):
            nuclei.coordinates = coordinate_list.item(i), coordinate_list.item(i+1), coordinate_list.item(i+2)
        basis_set = read_basis_set_file(self.basis_file, self.nuclei_list)
        with redirect_stdout(open(os.devnull, "w")):
            energy = self.energy_object.calculate_energy(self.nuclei_list, basis_set)
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
