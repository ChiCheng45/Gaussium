import numpy


class Matrix:

    def __init__(self, nuclei_array, file_reader_basis):
        self.nuclei_array = nuclei_array
        self.file_reader_basis = file_reader_basis

    def create_matrix(self, integrate):
        i = 0
        number_of_nuclei = len(self.nuclei_array)
        A = []
        while i < number_of_nuclei:
            j = 0
            row = []
            while j < number_of_nuclei:
                if j <= i:
                    A_ij = integrate.calculate_integral(self.nuclei_array, self.file_reader_basis, i, j)
                    row.append(A_ij)
                else:
                    row.append(0)
                j += 1
            A.append(row)
            i += 1
        A = numpy.matrix(A)
        A = A + A.T - numpy.diag(numpy.diag(A))
        return A
