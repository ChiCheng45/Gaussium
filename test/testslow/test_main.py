from unittest import TestCase
from numpy import testing
from src.main import start


class TestMain(TestCase):

    def test_total_rhf_energy_of_ethene_for_the_3_21_g_basis_set(self):
        total_energy = start('C2H4.mol', '3-21G.gbs', 'RHF', 1)
        testing.assert_approx_equal(total_energy, -77.600460844, 6)

    def test_total_rhf_energy_of_ethene_for_the_3_21_g_basis_set_with_symmetry(self):
        total_energy = start('C2H4.mol', '3-21G.gbs', 'RHF', 1, True)
        testing.assert_approx_equal(total_energy, -77.600460844, 6)

    def test_total_uhf_energy_of_oxygen_for_the_sto_3g_basis_set_with_symmetry(self):
        total_energy = start('O2.mol', 'STO-3G.gbs', 'UHF', 1, True)
        testing.assert_approx_equal(total_energy, -147.634028141, 6)

    def test_total_guhf_energy_of_oxygen_for_the_sto_3g_basis_set_with_symmetry(self):
        total_energy = start('O2.mol', 'STO-3G.gbs', 'GUHF', 1, True)
        testing.assert_approx_equal(total_energy, -147.634028141, 6)

    def test_total_mp2_energy_of_carbon_monoxide_for_the_sto_3g_basis_set(self):
        total_energy = start('CO.mol', 'STO-3G.gbs', 'MP2', 1)
        testing.assert_approx_equal(total_energy, -111.354512528, 6)

    def test_total_dft_s_energy_of_he_for_the_sto_3g_basis_set(self):
        total_energy = start('He.mol', 'STO-3G.gbs', ('DFT', 'S', ''), 1)
        testing.assert_approx_equal(total_energy, -2.657311972, 6)

    def test_total_dft_svwn3_energy_of_he_for_the_sto_3g_basis_set(self):
        total_energy = start('He.mol', 'STO-3G.gbs', ('DFT', 'S', 'VWN3'), 1)
        testing.assert_approx_equal(total_energy, -2.809598595, 6)
