from unittest import TestCase

from numpy import testing

from src.main import start


class TestMain(TestCase):

    def test_total_rhf_energy_of_heh_cation_for_the_sto_3g_basis_set(self):
        total_energy = start('HeH+.mol', 'STO-3G.gbs', 'RHF', 1)
        testing.assert_approx_equal(total_energy, -2.84183608212, 6)

    def test_total_rhf_energy_of_heh_cation_for_the_6_311_plus_g_star_star_basis_set(self):
        total_energy = start('HeH+.mol', '6-311+GPP.gbs', 'RHF', 1)
        testing.assert_approx_equal(total_energy, -2.92922773384, 6)
