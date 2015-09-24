from unittest import TestCase
from src.main import Nuclei
from random import random


class TestNuclei(TestCase):

    def setUp(self):
        self.a = random()
        self.b = random()
        self.c = random()
        self.d = random()
        self.e = random()
        self.nuclei = Nuclei(['HYDROGEN', self.a, self.b, self.c, self.d, self.e])

    def test_get_name(self):
        name = self.nuclei.get_name()
        self.assertEquals(name, 'HYDROGEN')

    def test_get_charge(self):
        charge = self.nuclei.get_charge()
        self.assertEquals(charge, self.a)

    def test_get_mass(self):
        mass = self.nuclei.get_mass()
        self.assertEquals(mass, self.b)

    def test_get_x(self):
        x = self.nuclei.get_x()
        self.assertEquals(x, self.c)

    def test_get_y(self):
        y = self.nuclei.get_y()
        self.assertEquals(y, self.d)

    def test_get_z(self):
        z = self.nuclei.get_z()
        self.assertEquals(z, self.e)
