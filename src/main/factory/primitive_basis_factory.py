from src.main.objects import PrimitiveBasis
from src.main.objects import Basis


def expand_basis_set(file_list, coordinates):
    basis_list = []
    for a in range(len(file_list)):
        if file_list[a][0] == 'S':
            primitive_basis_s_list = []
            for b in range(1, len(file_list[a])):
                primitive_basis = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 0, 0))
                primitive_basis_s_list.append(primitive_basis)
            basis_s = Basis(primitive_basis_s_list, coordinates, (0, 0, 0))
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
            basis_s = Basis(primitive_basis_s_list, coordinates, (0, 0, 0))
            basis_px = Basis(primitive_basis_px_list, coordinates, (1, 0, 0))
            basis_py = Basis(primitive_basis_py_list, coordinates, (0, 1, 0))
            basis_pz = Basis(primitive_basis_pz_list, coordinates, (0, 0, 1))
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
            basis_px = Basis(primitive_basis_px_list, coordinates, (1, 0, 0))
            basis_py = Basis(primitive_basis_py_list, coordinates, (0, 1, 0))
            basis_pz = Basis(primitive_basis_pz_list, coordinates, (0, 0, 1))
            basis_list.append(basis_px)
            basis_list.append(basis_py)
            basis_list.append(basis_pz)
        elif file_list[a][0] == 'D':
            primitive_basis_dxx_list = []
            primitive_basis_dyy_list = []
            primitive_basis_dzz_list = []
            primitive_basis_dxy_list = []
            primitive_basis_dxz_list = []
            primitive_basis_dyz_list = []
            for b in range(1, len(file_list[a])):
                primitive_basis_dxx = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (2, 0, 0))
                primitive_basis_dyy = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 2, 0))
                primitive_basis_dzz = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 0, 2))
                primitive_basis_dxy = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (1, 1, 0))
                primitive_basis_dxz = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (1, 0, 1))
                primitive_basis_dyz = PrimitiveBasis(file_list[a][b][0], file_list[a][b][1], coordinates, (0, 1, 1))
                primitive_basis_dxx_list.append(primitive_basis_dxx)
                primitive_basis_dyy_list.append(primitive_basis_dyy)
                primitive_basis_dzz_list.append(primitive_basis_dzz)
                primitive_basis_dxy_list.append(primitive_basis_dxy)
                primitive_basis_dxz_list.append(primitive_basis_dxz)
                primitive_basis_dyz_list.append(primitive_basis_dyz)
            basis_dxx = Basis(primitive_basis_dxx_list, coordinates, (2, 0, 0))
            basis_dyy = Basis(primitive_basis_dyy_list, coordinates, (0, 2, 0))
            basis_dzz = Basis(primitive_basis_dzz_list, coordinates, (0, 0, 2))
            basis_dxy = Basis(primitive_basis_dxy_list, coordinates, (1, 1, 0))
            basis_dxz = Basis(primitive_basis_dxz_list, coordinates, (1, 0, 1))
            basis_dyz = Basis(primitive_basis_dyz_list, coordinates, (0, 1, 1))
            basis_list.append(basis_dxx)
            basis_list.append(basis_dyy)
            basis_list.append(basis_dzz)
            basis_list.append(basis_dxy)
            basis_list.append(basis_dxz)
            basis_list.append(basis_dyz)
    return basis_list


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


def gaussian_product(gaussian1, gaussian2):
    d_1 = gaussian1.contraction
    a_1 = gaussian1.exponent
    l_1 = gaussian1.integral_exponents
    r_1 = gaussian1.coordinates
    d_2 = gaussian2.contraction
    a_2 = gaussian2.exponent
    l_2 = gaussian2.integral_exponents
    r_2 = gaussian2.coordinates

    d_3 = d_1 * d_2
    a_3 = a_1 + a_2
    l_3 = (l_1[0] + l_2[0], l_1[1] + l_2[1], l_1[2] + l_2[2])
    r_3 = (a_1*r_1[0] + a_2*r_2[0]) / a_3, (a_1*r_1[1] + a_2*r_2[1]) / a_3, (a_1*r_1[2] + a_2*r_2[2]) / a_3
    return PrimitiveBasis(d_3, a_3, r_3, l_3)
