import numpy


class Matrix:

    def __init__(self, nuclei_array):
        self.nuclei_array = nuclei_array

    def create_matrix(self, element):
        number_of_nuclei = len(self.nuclei_array)
        matrix = []
        for i in range(0, number_of_nuclei):
            j = 0
            row = []
            for j in range(0, number_of_nuclei):
                if j <= i:
                    matrix_ij = element.calculate(i, j)
                    row.append(matrix_ij)
                else:
                    row.append(0)
            matrix.append(row)
        matrix = numpy.matrix(matrix)
        matrix = matrix + matrix.T - numpy.diag(numpy.diag(matrix))
        return matrix
