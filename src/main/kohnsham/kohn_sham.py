from src.main.hartreefock import Restricted
from src.main.hartreefock import RestrictedSCF
from src.main.kohnsham import slater_exchange
from src.main.kohnsham import chachiyo_correlation
from src.main.kohnsham import ExchangeCorrelation
from src.main.kohnsham import RestrictedKohnShamHamiltonian


class RestrictedKohnSham(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes, exchange, correlation):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

        if exchange == 'SLATER':
            exchange = slater_exchange
        else:
            exchange = lambda x: 0

        if correlation == 'CHACHIYO':
            correlation = chachiyo_correlation
        else:
            correlation = lambda x: 0

        self.scf_method = RestrictedSCF(self.linear_algebra, self.electrons, self.orbital_overlap,
        RestrictedKohnShamHamiltonian(self.core_hamiltonian, self.repulsion,
        ExchangeCorrelation(basis_set_array, exchange, correlation)))
        print('\n\nBEGIN RESTRICTED KOHN SHAM\n')
