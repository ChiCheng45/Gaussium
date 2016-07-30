from src.main.hartreefock import Restricted
from src.main.hartreefock import RestrictedSCF
from src.main.kohnsham import SlaterExchange
from src.main.kohnsham import RestrictedKohnShamHamiltonian


class RestrictedKohnSham(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes, exchange_functional):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

        if exchange_functional == 'Slater':
            exchange_functional = SlaterExchange(basis_set_array)

        self.scf_method = RestrictedSCF(self.linear_algebra, self.electrons, self.orbital_overlap,
        RestrictedKohnShamHamiltonian(self.core_hamiltonian, self.repulsion, exchange_functional))
        print('\n\nBEGIN RESTRICTED KOHN SHAM\n')
