from src.main.hartreefock import Restricted
from src.main.hartreefock import RestrictedSCF
from src.main.kohnsham.exchange import ExchangePotential
from src.main.kohnsham.exchange import SlaterExchange
from src.main.kohnsham.correlation import CorrelationPotential
from src.main.kohnsham.correlation import VoskoWilkNusair
from src.main.kohnsham import ExchangeCorrelation
from src.main.kohnsham import RestrictedKohnShamHamiltonian


class RestrictedKohnSham(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes, exchange, correlation):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

        if exchange == 'S':
            exchange = SlaterExchange(alpha=1.0)
        elif exchange == 'XA':
            exchange = SlaterExchange(alpha=0.7)
        else:
            exchange = ExchangePotential()  # returns a potential of 0.0

        if correlation == 'VWN3':
            correlation = VoskoWilkNusair(a=0.0621814, x_0=-0.409286, b=13.0720, c=42.7198)
        elif correlation == 'VWN5':
            correlation = VoskoWilkNusair(a=0.0621814, x_0=-0.10498, b=3.72744, c=12.9352)
        else:
            correlation = CorrelationPotential()  # returns a potential of 0.0

        self.scf_method = RestrictedSCF(
            self.linear_algebra, self.electrons, self.orbital_overlap, RestrictedKohnShamHamiltonian(
                self.core_hamiltonian, self.repulsion, ExchangeCorrelation(basis_set_array, exchange, correlation)
            )
        )

        print('\n\nBEGIN RESTRICTED KOHN SHAM\n')
