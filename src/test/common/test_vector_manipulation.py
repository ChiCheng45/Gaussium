from unittest import TestCase
from numpy import testing
from math import pi
from src.main.common import Vector


class TestVectorQuaternions(TestCase):

    def test_create_quaternion_returns_quaternion_for_0_degree_rotation_around_x_axis(self):
        axis_of_rotation = (1.0, 0.0, 0.0)
        theta = 0
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (1.0, 0.0, 0.0, 0.0), 6)

    def test_create_quaternion_returns_quaternion_for_90_degree_rotation_around_x_axis(self):
        axis_of_rotation = (1.0, 0.0, 0.0)
        theta = pi / 2
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (0.707107, 0.707107, 0.0, 0.0), 6)

    def test_create_quaternion_returns_quaternion_for_180_degree_rotation_around_x_axis(self):
        axis_of_rotation = (1.0, 0.0, 0.0)
        theta = pi
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (0.0, 1.0, 0.0, 0.0), 6)

    def test_create_quaternion_returns_quaternion_for_90_degree_rotation_around_y_axis(self):
        axis_of_rotation = (0.0, 1.0, 0.0)
        theta = pi / 2
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (0.707107, 0.0, 0.707107, 0.0), 6)

    def test_create_quaternion_returns_quaternion_for_180_degree_rotation_around_y_axis(self):
        axis_of_rotation = (0.0, 1.0, 0.0)
        theta = pi
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (0.0, 0.0, 1.0, 0.0), 6)

    def test_create_quaternion_returns_quaternion_for_45_degree_rotation_around_x_axis(self):
        axis_of_rotation = (1.0, 0.0, 0.0)
        theta = pi / 4
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        testing.assert_array_almost_equal(quaternion, (0.923880, 0.382683, 0.0, 0.0), 6)

    def test_rotation_of_unit_y_vector_with_axis_of_rotation_y_by_90_degrees(self):
        theta = pi / 2
        axis_of_rotation = (0, 1, 0)
        x_vector = (1, 0, 0)
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        rotated_vector = Vector.quaternion_rotation(quaternion, x_vector)

        testing.assert_array_almost_equal(rotated_vector, (0, 0, -1), 6)

    def test_rotation_of_unit_z_vector_with_axis_of_rotation_x_by_45_degrees(self):
        theta = pi / 4
        axis_of_rotation = (1, 0, 0)
        z_vector = (0, 0, -1)
        quaternion = Vector.create_quaternion(axis_of_rotation, theta)
        rotated_vector = Vector.quaternion_rotation(quaternion, z_vector)
        testing.assert_array_almost_equal(rotated_vector, (0, 0.707107, -0.707107), 6)
