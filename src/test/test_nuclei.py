from unittest import TestCase
from src.main import Nuclei
from random import random


class TestNuclei(TestCase):

    def setUp(self):
        self.nuclei = Nuclei(['HYDROGEN', 1, 3, 5, 7, 11])

    def test_get_name(self):
        name = self.nuclei.get_name()
        self.assertEquals(name, 'HYDROGEN')

    def test_get_charge(self):
        charge = self.nuclei.get_charge()
        self.assertEquals(charge, 1)

    def test_get_mass(self):
        mass = self.nuclei.get_mass()
        self.assertEquals(mass, 3)

    def test_get_x(self):
        x = self.nuclei.get_x()
        self.assertEquals(x, 5)

    def test_get_y(self):
        y = self.nuclei.get_y()
        self.assertEquals(y, 7)

    def test_get_z(self):
        z = self.nuclei.get_z()
        self.assertEquals(z, 11)
