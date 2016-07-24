class CoupledCluster:

    def __init__(self, hartree_fock):
        self.hartree_fock = hartree_fock

    def singles_doubles(self):
        electron_energy, orb_energies, orbital_coefficients = self.hartree_fock.begin_scf()
