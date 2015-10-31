from src.main.objects import PrimitiveBasis

class PrimitiveBasisFactory:

    def expand_basis(self, function_type, coefficients_array):
        primitive_array = []
        if function_type == 'S':
            primitive_basis = PrimitiveBasis(function_type, coefficients_array[0], coefficients_array[1], [0, 0, 0])
            primitive_array.append(primitive_basis)
            return primitive_array
        elif function_type == 'L':
            return primitive_array

    def del_operator(self, primitive_gaussian):
        array = []
        orbital = primitive_gaussian.get_orbital_type()
        l = primitive_gaussian.get_integral_exponents()[0]
        m = primitive_gaussian.get_integral_exponents()[1]
        n = primitive_gaussian.get_integral_exponents()[2]
        a = primitive_gaussian.get_exponent()
        c = primitive_gaussian.get_contraction()

        contraction = a * ((2 * (l + m + n)) + 3)
        primitive_basis = PrimitiveBasis(orbital, contraction, a, [l, m, n])
        array.append(primitive_basis)

        contraction = - 2 * a**2
        primitive_basis = PrimitiveBasis(orbital, contraction, a, [l + 2, m, n])
        array.append(primitive_basis)
        primitive_basis = PrimitiveBasis(orbital, contraction, a, [l, m + 2, n])
        array.append(primitive_basis)
        primitive_basis = PrimitiveBasis(orbital, contraction, a, [l, m, n + 2])
        array.append(primitive_basis)

        if l >= 2:
            contraction = - (1/2) * (l * (l - 1))
            primitive_basis = PrimitiveBasis(orbital, contraction, a, [l - 2, m, n])
            array.append(primitive_basis)
        if m >= 2:
            contraction = - (1/2) * (m * (m - 1))
            primitive_basis = PrimitiveBasis(orbital, contraction, a, [l, m - 2, n])
            array.append(primitive_basis)
        if n >= 2:
            contraction = - (1/2) * (n * (n - 1))
            primitive_basis = PrimitiveBasis(orbital, contraction, a, [l, m, n - 2])
            array.append(primitive_basis)

        return array
