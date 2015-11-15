from src.main.objects import PrimitiveBasis
from src.main.objects import Basis


class PrimitiveBasisFactory:

    @staticmethod
    def expand_basis(file_list, coordinates, nuclei):
        basis_list = []
        for a in range(len(file_list)):
            if file_list[a][0] == 'S':
                primitive_basis_s_list = []
                for b in range(1, len(file_list[a])):
                    primitive_basis = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 0, 0))
                    primitive_basis_s_list.append(primitive_basis)
                basis_s = Basis(primitive_basis_s_list)
                basis_list.append(basis_s)
            elif file_list[a][0] == 'L':
                primitive_basis_s_list = []
                primitive_basis_px_list = []
                primitive_basis_py_list = []
                primitive_basis_pz_list = []
                for b in range(1, len(file_list[a])):
                    primitive_basis_s = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 0, 0))
                    primitive_basis_px = PrimitiveBasis(file_list[a][b][2], file_list[a][b][1], coordinates, (1, 0, 0))
                    primitive_basis_py = PrimitiveBasis(file_list[a][b][2], file_list[a][b][1], coordinates, (0, 1, 0))
                    primitive_basis_pz = PrimitiveBasis(file_list[a][b][2], file_list[a][b][1], coordinates, (0, 0, 1))
                    primitive_basis_s_list.append(primitive_basis_s)
                    primitive_basis_px_list.append(primitive_basis_px)
                    primitive_basis_py_list.append(primitive_basis_py)
                    primitive_basis_pz_list.append(primitive_basis_pz)
                basis_s = Basis(primitive_basis_s_list)
                basis_px = Basis(primitive_basis_px_list)
                basis_py = Basis(primitive_basis_py_list)
                basis_pz = Basis(primitive_basis_pz_list)
                basis_list.append(basis_s)
                basis_list.append(basis_px)
                basis_list.append(basis_py)
                basis_list.append(basis_pz)
            elif file_list[a][0] == 'P':
                primitive_basis_px_list = []
                primitive_basis_py_list = []
                primitive_basis_pz_list = []
                for b in range(1, len(file_list[a])):
                    primitive_basis_px = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (1, 0, 0))
                    primitive_basis_py = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 1, 0))
                    primitive_basis_pz = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 0, 1))
                    primitive_basis_px_list.append(primitive_basis_px)
                    primitive_basis_py_list.append(primitive_basis_py)
                    primitive_basis_pz_list.append(primitive_basis_pz)
                basis_px = Basis(primitive_basis_px_list)
                basis_py = Basis(primitive_basis_py_list)
                basis_pz = Basis(primitive_basis_pz_list)
                basis_list.append(basis_px)
                basis_list.append(basis_py)
                basis_list.append(basis_pz)
        return basis_list

    @staticmethod
    def del_operator(primitive_gaussian):
        primitive_array = []
        l = primitive_gaussian.integral_exponents[0]
        m = primitive_gaussian.integral_exponents[1]
        n = primitive_gaussian.integral_exponents[2]
        a = primitive_gaussian.exponent
        coordinates = primitive_gaussian.coordinates

        contraction = a * ((2 * (l + m + n)) + 3)
        primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l, m, n))
        primitive_array.append(primitive_basis)

        contraction = - 2 * a**2
        primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l + 2, m, n))
        primitive_array.append(primitive_basis)
        primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l, m + 2, n))
        primitive_array.append(primitive_basis)
        primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l, m, n + 2))
        primitive_array.append(primitive_basis)

        if l >= 2:
            contraction = - (1/2) * (l * (l - 1))
            primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l - 2, m, n))
            primitive_array.append(primitive_basis)
        if m >= 2:
            contraction = - (1/2) * (m * (m - 1))
            primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l, m - 2, n))
            primitive_array.append(primitive_basis)
        if n >= 2:
            contraction = - (1/2) * (n * (n - 1))
            primitive_basis = PrimitiveBasis(contraction, a, coordinates, (l, m, n - 2))
            primitive_array.append(primitive_basis)

        return primitive_array
