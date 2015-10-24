from src.main.objects import PrimitiveBasis

class PrimitiveBasisFactory:

    def expand_basis(self, function_type, coefficients_array):
        primitive_array = []
        if function_type == 'S':
            primitive_basis = PrimitiveBasis(function_type, coefficients_array[0], coefficients_array[1], [0, 0, 0])
            primitive_array.append(primitive_basis)
            return primitive_array
        elif function_type == 'L':
            return []
