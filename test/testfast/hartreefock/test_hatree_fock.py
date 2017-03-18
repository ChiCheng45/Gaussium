from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.common import Symmetry
from src.hartreefock import RestrictedHF


class TestHartreeFock(TestCase):

    def setUp(self):
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000, 0.0000, 1.4632))
        helium_1 = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(0.0000, 0.0000, 0.0000))

        hydrogen_s_1 = MagicMock(
            contraction=0.15432897, exponent=3.42525091, coordinates=(0.0000, 0.0000, 1.4632),
            integral_exponents=(0, 0, 0), normalisation=1.794441832218435
        )
        hydrogen_s_2 = MagicMock(
            contraction=0.53532814, exponent=0.62391373, coordinates=(0.0000, 0.0000, 1.4632),
            integral_exponents=(0, 0, 0), normalisation=0.5003264923314032
        )
        hydrogen_s_3 = MagicMock(
            contraction=0.44463454, exponent=0.16885540, coordinates=(0.0000, 0.0000, 1.4632),
            integral_exponents=(0, 0, 0), normalisation=0.18773545851092535
        )

        helium_s_1 = MagicMock(
            contraction=0.15432897, exponent=9.75393461, coordinates=(0.0000, 0.0000, 0.0000),
            integral_exponents=(0, 0, 0), normalisation=3.9336432656254527
        )
        helium_s_2 = MagicMock(
            contraction=0.53532814, exponent=1.77669115, coordinates=(0.0000, 0.0000, 0.0000),
            integral_exponents=(0, 0, 0), normalisation=1.0967787981767012
        )
        helium_s_3 = MagicMock(
            contraction=0.44463454, exponent=0.48084429, coordinates=(0.0000, 0.0000, 0.0000),
            integral_exponents=(0, 0, 0), normalisation=0.41154131374122654
        )

        hydrogen_basis = MagicMock(
            primitive_gaussian_array=[hydrogen_s_1, hydrogen_s_2, hydrogen_s_3], coordinates=(0.0000, 0.0000, 1.4632),
            integral_exponents=(0, 0, 0)
        )
        helium_basis = MagicMock(
            primitive_gaussian_array=[helium_s_1, helium_s_2, helium_s_3], coordinates=(0.0000, 0.0000, 0.0000),
            integral_exponents=(0, 0, 0)
        )

        point_group = MagicMock(
            rotation_symmetry=[], reflection_symmetry=[], improper_rotation=[], inversion_symmetry=[], label='C_{1}'
        )

        nuclei_array_heh = [hydrogen_1, helium_1]
        basis_set = [hydrogen_basis, helium_basis]
        symmetry_object = Symmetry(point_group, basis_set)

        self.restricted_hartree_fock = RestrictedHF(nuclei_array_heh, basis_set, 2, symmetry_object, 1)

    def test_electronic_energy_of_heh_cation_for_the_sto_3g_basis_set(self):
        electron_energy = self.restricted_hartree_fock.begin_scf()[0]
        testing.assert_approx_equal(electron_energy, -4.227529, 6)
