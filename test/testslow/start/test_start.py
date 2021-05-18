import os
from unittest import TestCase
from numpy import testing
from gaussium.main import start


class TestStart(TestCase):

    def test_total_rhf_energy_of_heh_cation_for_the_sto_3g_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, 'RHF', 1)
        testing.assert_approx_equal(total_energy, -2.84183608212, 6)

    def test_total_rhf_energy_of_heh_cation_for_the_6_311_plus_g_star_star_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'HeH+.mol')
        basis_file = os.path.join(os.path.dirname(__file__), '6-311+GPP.gbs')
        total_energy = start(str(mol_file), str(basis_file), 'RHF', 1)
        testing.assert_approx_equal(total_energy, -2.92922773384, 6)

    def test_total_rhf_energy_of_ethene_for_the_3_21_g_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'C2H4.mol')
        basis_file = os.path.join(os.path.dirname(__file__), '3-21G.gbs')
        total_energy = start(mol_file, basis_file, 'RHF', 1)
        testing.assert_approx_equal(total_energy, -77.600460844, 6)

    def test_total_rhf_energy_of_ethene_for_the_3_21_g_basis_set_with_symmetry(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'C2H4.mol')
        basis_file = os.path.join(os.path.dirname(__file__), '3-21G.gbs')
        total_energy = start(mol_file, basis_file, 'RHF', 1, True)
        testing.assert_approx_equal(total_energy, -77.600460844, 6)

    def test_total_uhf_energy_of_oxygen_for_the_sto_3g_basis_set_with_symmetry(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'O2.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, 'UHF', 1, True)
        testing.assert_approx_equal(total_energy, -147.634028141, 6)

    def test_total_guhf_energy_of_oxygen_for_the_sto_3g_basis_set_with_symmetry(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'O2.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, 'GUHF', 1, True)
        testing.assert_approx_equal(total_energy, -147.634028141, 6)

    def test_total_mp2_energy_of_carbon_monoxide_for_the_sto_3g_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'CO.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, 'MP2', 1)
        testing.assert_approx_equal(total_energy, -111.354512528, 6)

    def test_total_dft_s_energy_of_he_for_the_sto_3g_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'He.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, ('DFT', 'S', ''), 1)
        testing.assert_approx_equal(total_energy, -2.657311972, 6)

    def test_total_dft_svwn3_energy_of_he_for_the_sto_3g_basis_set(self):
        mol_file = os.path.join(os.path.dirname(__file__), 'He.mol')
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        total_energy = start(mol_file, basis_file, ('DFT', 'S', 'VWN3'), 1)
        testing.assert_approx_equal(total_energy, -2.809598595, 6)
