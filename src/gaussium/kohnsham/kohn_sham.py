from gaussium.kohnsham import ExchangeCorrelation
from gaussium.kohnsham import RestrictedKohnShamHamiltonian
from gaussium.kohnsham.exchange import Exchange
from gaussium.kohnsham.exchange import SlaterExchange
from gaussium.hartreefock import Restricted
from gaussium.kohnsham.kohn_sham_scf import KSRestrictedSCF
from gaussium.kohnsham.correlation import Correlation
from gaussium.kohnsham.correlation import VoskoWilkNusair


class RestrictedKohnSham(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes, exchange, correlation):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

        if exchange == 'S':
            exchange = SlaterExchange(alpha=2/3)
        elif exchange == 'XA':
            exchange = SlaterExchange(alpha=0.7)
        else:
            exchange = Exchange()  # returns a potential of 0.0

        if correlation == 'VWN3':
            correlation = VoskoWilkNusair(a=0.0621814, x_0=-0.409286, b=13.0720, c=42.7198)
        elif correlation == 'VWN5':
            correlation = VoskoWilkNusair(a=0.0621814, x_0=-0.10498, b=3.72744, c=12.9352)
        else:
            correlation = Correlation()  # returns a potential of 0.0

        xc = ExchangeCorrelation(basis_set_array, exchange, correlation)
        self.scf_method = KSRestrictedSCF(
            self.linear_algebra, self.electrons, self.orbital_overlap, xc,
            RestrictedKohnShamHamiltonian(
                self.core_hamiltonian, self.repulsion, xc
            )
        )

        print('\n\nBEGIN RESTRICTED KOHN SHAM\n')
