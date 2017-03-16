from src.objects import ImproperRotationSymmetry
from src.objects import InversionSymmetry
from src.objects import ReflectionSymmetry
from src.objects import RotationSymmetry


class PointGroup:

    def __init__(self, rotation_symmetry, reflection_symmetry, improper_rotation, inversion_symmetry, label):
        self.rotation_symmetry = rotation_symmetry
        self.reflection_symmetry = reflection_symmetry
        self.improper_rotation = improper_rotation
        self.inversion_symmetry = inversion_symmetry
        self.label = label


class Oh(PointGroup):

    def __init__(self):
        rotation_symmetry = [
            RotationSymmetry(4, (0.0, 0.0, 1.0)),
            RotationSymmetry(4, (0.0, 1.0, 0.0)),
            RotationSymmetry(4, (1.0, 0.0, 0.0)),
            RotationSymmetry(2, (0.0, 0.707107, 0.707107)),
            RotationSymmetry(2, (0.0, -0.707107, 0.707107)),
            RotationSymmetry(2, (0.707107, 0.0, 0.707107)),
            RotationSymmetry(2, (-0.707107, 0.0, 0.707107)),
            RotationSymmetry(2, (0.707107, 0.707107, 0.0)),
            RotationSymmetry(2, (-0.707107, 0.707107, 0.0)),
            RotationSymmetry(3, (0.557735, 0.557735, 0.557735)),
            RotationSymmetry(3, (-0.557735, 0.557735, 0.557735)),
            RotationSymmetry(3, (0.557735, -0.557735, 0.557735)),
            RotationSymmetry(3, (-0.557735, -0.557735, 0.557735)),
        ]

        reflection_symmetry = [
            ReflectionSymmetry((0.0, 0.0, 1.0)),
            ReflectionSymmetry((0.0, 1.0, 0.0)),
            ReflectionSymmetry((1.0, 0.0, 0.0)),
            ReflectionSymmetry((0.0, 0.707107, 0.707107)),
            ReflectionSymmetry((0.0, -0.707107, 0.707107)),
            ReflectionSymmetry((0.707107, 0.0, 0.707107)),
            ReflectionSymmetry((-0.707107, 0.0, 0.707107)),
            ReflectionSymmetry((0.707107, 0.707107, 0.0)),
            ReflectionSymmetry((-0.707107, 0.707107, 0.0))
        ]

        improper_rotation = [
            ImproperRotationSymmetry(6, (0.557735, 0.557735, 0.557735)),
            ImproperRotationSymmetry(6, (-0.557735, 0.557735, 0.557735)),
            ImproperRotationSymmetry(6, (0.557735, -0.557735, 0.557735)),
            ImproperRotationSymmetry(6, (-0.557735, -0.557735, 0.557735)),
            ImproperRotationSymmetry(4, (0.0, 0.0, 1.0)),
            ImproperRotationSymmetry(4, (0.0, 1.0, 0.0)),
            ImproperRotationSymmetry(4, (1.0, 0.0, 0.0))
        ]

        inversion_symmetry = [
            InversionSymmetry()
        ]

        label = 'O_{h}'

        super().__init__(rotation_symmetry, reflection_symmetry, improper_rotation, inversion_symmetry, label)


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
            ReflectionSymmetry((0.0, 0.0, 1.0)),
            ReflectionSymmetry((0.0, 1.0, 0.0)),
            ReflectionSymmetry((1.0, 0.0, 0.0)),
            ReflectionSymmetry((0.707107, 0.707107, 0.0)),
            ReflectionSymmetry((-0.707107, 0.707107, 0.0))
        ]

        improper_rotation = [
            ImproperRotationSymmetry(4, (0.0, 0.0, 1.0))
        ]

        inversion_symmetry = [
            InversionSymmetry()
        ]

        label = 'D_{4h}'

        super().__init__(rotation_symmetry, reflection_symmetry, improper_rotation, inversion_symmetry, label)


class C4v(PointGroup):

    def __init__(self):
        rotation_symmetry = [
            RotationSymmetry(4, (0.0, 0.0, 1.0))
        ]

        reflection_symmetry = [
            ReflectionSymmetry((0.0, 1.0, 0.0)),
            ReflectionSymmetry((1.0, 0.0, 0.0)),
            ReflectionSymmetry((0.707107, 0.707107, 0.0)),
            ReflectionSymmetry((-0.707107, 0.707107, 0.0))
        ]

        improper_rotation = []

        inversion_symmetry = []

        label = 'C_{4v}'

        super().__init__(rotation_symmetry, reflection_symmetry, improper_rotation, inversion_symmetry, label)
