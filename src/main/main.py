from src.main.hartreefock import HartreeFock

if __name__ == "__main__":
    HartreeFock.begin('O2.mol', 'STO-3G.gbs', 'UHF')
    # HartreeFock.begin('O2.mol', 'STO-3G.gbs', 'RHF')
