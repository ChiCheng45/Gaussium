from gaussium import start


def run():
    # example calculations here we use the backslashes for the path
    # to the molecule and basis set files since we are on windows
    # change to forward slash for unix

    # start('molfiles\\HeH+.mol', 'basissets\\STO-3G.gbs', 'RHF', 4)  # -2.84183608212 a.u. (energy from reference calculation)
    # start('molfiles\\HeH+.mol', 'basissets\\6-311+GPP.gbs', 'RHF', 4)  # -2.92922773384 a.u.
    # start('molfiles\\C2H4.mol', 'basissets\\3-21G.gbs', 'RHF', 4)  # -77.600460844 a.u. 30.747198048700866s
    # start('molfiles\\O2.mol', 'basissets\\STO-3G.gbs', 'UHF', 4, symmetry=True)  # -147.634028141 a.u.
    # start('molfiles\\O2.mol', 'basissets\\STO-3G.gbs', 'GUHF', 4)  # -147.634028141 a.u.
    # start('molfiles\\CO.mol', 'basissets\\STO-3G.gbs', 'MP2', 4)  # -111.354512528 a.u.
    # start('molfiles\\H2O.mol', 'basissets\\STO-3G.gbs', 'RHF', 4, symmetry=True)
    # start('molfiles\\C2H4.mol', 'basissets\\3-21G.gbs', 'RHF', 4, symmetry=True)  # -77.600460844 a.u. 19.0269839632222s
    # start('molfiles\\H2O.mol', 'basissets\\STO-3G.gbs', 'CIS', 4)  # 0.2872554996 a.u. 0.3564617587 a.u.
    # start('molfiles\\He.mol', 'basissets\\3-21G.gbs', 'RHF', 4)  # -2.83567987364 a.u.
    # start('molfiles\\He.mol', 'basissets\\6-311G.gbs', 'RHF', 4)  # -2.85989542457 a.u.
    # start('molfiles\\He.mol', 'basissets\\cc-pVDZ.gbs', 'RHF', 4)  # -2.85516047724192 a.u.

    # start('molfiles\\H2O.mol', 'basissets\\STO-3G.gbs', 'CCSD', 4)  # -0.0706800939192 a.u.
    # start('molfiles\\CH4.mol', 'basissets\\STO-3G.gbs', 'CCSD', 4)  # -0.078469894846414 a.u.
    # start('molfiles\\H2O.mol', 'basissets\\STO-3G.gbs', 'CCSD(T)', 4)  # -9.98772699528e-05 a.u.

    # geometry optimization
    # start('molfiles\\H2O.mol', 'basissets\\STO-3G.gbs', 'RHF', 4, geometry_optimization='NelderMead')  # -74.96588377357489 a.u.

    # only worth doing DFT calculations on atoms at the moment
    # start('molfiles\\He.mol', 'basissets\\STO-3G.gbs', ('DFT', 'S', ''), 4)  # -2.65731197167 a.u.
    # start('molfiles\\He.mol', 'basissets\\3-21G.gbs', ('DFT', 'S', ''), 4)  # -2.69341499501 a.u.
    # start('molfiles\\He.mol', 'basissets\\STO-6G.gbs', ('DFT', 'S', ''), 4)  # -2.69600757420 a.u.
    # start('molfiles\\Li-.mol', 'basissets\\STO-3G.gbs', ('DFT', 'S', ''), 4) # -6.94823326080 a.u.
    # start('molfiles\\He.mol', 'basissets\\STO-3G.gbs', ('DFT', 'XA', ''), 4)  # -2.70257401931 a.u.
    # start('molfiles\\He.mol', 'basissets\\STO-3G.gbs', ('DFT', 'S', 'VWN3'), 4)  # -2.80959859438 a.u.
    # start('molfiles\\Li-.mol', 'basissets\\STO-3G.gbs', ('DFT', 'S', 'VWN3'), 4)  # -7.22285707872 a.u.
    # start('molfiles\\He.mol', 'basissets\\3-21G.gbs', ('DFT', 'S', 'VWN3'), 4)  # -2.84354346745 a.u.
    # start('molfiles\\He.mol', 'basissets\\cc-pVDZ.gbs', ('DFT', 'S', 'VWN3'), 4)  # -2.86415519051 a.u.
    start('molfiles\\He.mol', 'basissets\\3-21G.gbs', ('DFT', 'S', 'VWN5'), 4)  # -2.80601675458 a.u.


if __name__ == "__main__":
    run()
