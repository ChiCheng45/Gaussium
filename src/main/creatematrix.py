import numpy


class Matrix:

    def __init__(self, nuclei_array):
        self.nuclei_array = nuclei_array

    def create_matrix(self, calc_element):
        i = 0
        number_of_nuclei = len(self.nuclei_array)
        matrix = []
        while i < number_of_nuclei:
            j = 0
            row = []
            while j < number_of_nuclei:
                if j <= i:
                    matrix_ij = calc_element.calculate(i, j)
                    row.append(matrix_ij)
                else:
                    row.append(0)
                j += 1
            matrix.append(row)
            i += 1
        matrix = numpy.matrix(matrix)
        matrix = matrix + matrix.T - numpy.diag(numpy.diag(matrix))
        return matrix
