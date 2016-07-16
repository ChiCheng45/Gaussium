from src.main.objects import RotationSymmetry
from src.main.objects import ReflectionSymmetry
from src.main.objects import InversionSymmetry


class PointGroup:

    def __init__(self, rotation_symmetry, reflection_symmetry, inversion_symmetry, label):
        self.rotation_symmetry = rotation_symmetry
        self.reflection_symmetry = reflection_symmetry
        self.inversion_symmetry = inversion_symmetry
        self.label = label


class D4h(PointGroup):

    def __init__(self):
        rotation_symmetry = [
            RotationSymmetry(4, (0.0, 0.0, 1.0)),
            RotationSymmetry(2, (0.0, 1.0, 0.0)),
            RotationSymmetry(2, (1.0, 0.0, 0.0)),
            RotationSymmetry(2, (0.707107, 0.707107, 0.0)),
            RotationSymmetry(2, (-0.707107, 0.707107, 0.0))
        ]

        reflection_symmetry = [
            ReflectionSymmetry((1.0, 0.0, 0.0)),
            ReflectionSymmetry((0.0, 1.0, 0.0)),
            ReflectionSymmetry((1.0, 0.0, 0.0)),
            ReflectionSymmetry((0.707107, 0.707107, 0.0)),
            ReflectionSymmetry((-0.707107, 0.707107, 0.0))
        ]

        inversion_symmetry = [
            InversionSymmetry()
        ]

        label = 'D_{4h}'

        super().__init__(rotation_symmetry, reflection_symmetry, inversion_symmetry, label)
