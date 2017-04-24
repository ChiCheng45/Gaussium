## Introduction
A basic quantum chemical program written in Python 3 using the numpy and scipy libraries.

Currently this program fully supports RHF, UHF, CIS, TDHF, DFT, CCSD and CCSD(T). My next plans are to implementing more CI based methods and work on speeding up the CCSD iterations. 

I used Attlia Szabo and Neil S. Ostlunds _Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory_ and David B. Cooks _Handbook of Computational Quantum Chemistry_ as my main references for the theories and methods behind the electronic structure calculations. The developers resources at http://www.psicode.org/developers.php were also invaluable to the success of project and had a number of excellent tutorials and programming examples.

## Instructions
* To run this program add the desired .mol and .gbs files to the `molfiles` and `basisset` directories.
* Next edit the `src/main/main.py` `menu()` function so that the desired calculations are made for example,
```python
def menu():
    start('H2O.mol', 'STO-3G.gbs', 'CCSD', 4)
```
for DFT calculation the functional are given inputted using a tuple for SVWN3,
```python
def menu():
    start('He.mol', 'STO-3G.gbs', ('DFT', 'S', 'VWN3'), 4)
```
the start function contains more options such as the number of processes used during the multiprocessing sections of the code and whether symmetry is turned on for faster integral calculations. See `start()` for more details,
```python
def start(mol, basis, method, processes, symmetry=False)
```
* Now run the main.py, for example on a Windows machine,
```
C:\Anaconda3\python.exe C:\Users\username\PycharmProjects\Quantum_Chemistry\src\main.py
```

## Supported Features and Methods
* Restricted Hartree-Fock
* Unrestricted Hartree-Fock
* Density Functional Theory - SVWN and S_X
* Møller–Plesset Second Order
* Coupled Cluster Singles and Doubles
* Coupled Cluster Singles and Doubles with Perturbation Triples
* Configuration Interaction Singles
* Time-Dependent Hartree-Fock
* DIIS for SCF Calculations
* Multiprocessing during ERI evaluations
* Reduced ERI evaluations with Symmetry
* Obara-Saika recursion relation for ERI
* Nelder-Mead method for geometry optimization

## Example Calculations and Energy Comparisons

### Restricted Hartree-Fock
##### HeH<sup>+</sup> STO-3G 
Completed my first calculation of HeH+ with a bond-length of 1.4632 a<sub>0</sub> using the STO-3G basis set. Comparison with Spartan Student Edition v5,
```
SCF model:
 A restricted Hartree-Fock SCF calculation will be
 performed using Pulay DIIS + Geometric Direct Minimization

 SCF total energy:      -2.8418365 hartrees

  Reason for exit: Successful completion 
  Quantum Calculation CPU Time :          .39
  Quantum Calculation Wall Time:         2.85
```
```
ORBITAL COEFFICIENTS
[[-0.87660574  0.79774813]
 [-0.20247895 -1.16783645]]

TOTAL ENERGY: -2.841836483 a.u.


Time Taken: 0.10287352007652968s
```

##### HeH<sup>+</sup> 6-311+G**
```
    SPARTAN STUDENT Quantum Mechanics Program:  (PC/x86)   Release  5.0.0v4

Job type: Single point.
Method: RHF
Basis set: 6-311+G**
Number of shells: 8
Number of basis functions: 12
Charge :     +1 
Multiplicity: 1

SCF model:
 A restricted Hartree-Fock SCF calculation will be
 performed using Pulay DIIS + Geometric Direct Minimization

 SCF total energy:      -2.9292278 hartrees
```
```
TOTAL NUCLEAR REPULSION ENERGY: 1.36687066232 a.u.
TOTAL ENERGY: -2.92922773418 a.u.

*********************************************************************************************************

Time Taken: 10.986138812057659s
```

##### C<sub>2</sub>H<sub>4</sub> 3-21G
```
Job type: Single point.
Method: RHF
Basis set: 3-21G(*)
Number of shells: 14
Number of basis functions: 26
Multiplicity: 1

SCF model:
 A restricted Hartree-Fock SCF calculation will be
 performed using Pulay DIIS + Geometric Direct Minimization

 SCF total energy:     -77.6009882 hartrees
```
```
TOTAL NUCLEAR REPULSION ENERGY: 33.4424184132 a.u.
TOTAL ENERGY: -77.6004608443 a.u.

*********************************************************************************************************

Time Taken: 72.70841306778921s
```

### Unrestricted Hartree-Fock
##### O<sub>2</sub> STO-3G
Completed UHF calculation of O<sub>2</sub> with a bond-length of 2.28541 a<sub>0</sub> using the STO-3G basis set. Comparison with Psi4,
```
    Alpha Virtual:                                                        

       3B1u    0.687444  

    Beta Occupied:                                                        

       1B1u  -20.409351     1Ag   -20.409200     2Ag    -1.488240  
       2B1u   -0.897684     3Ag    -0.554936     1B2u   -0.443771  
       1B3u   -0.443771  

    Beta Virtual:                                                         

       1B2g    0.270958     1B3g    0.270958     3B1u    0.777260  

    Final Occupation by Irrep:
             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u 
    DOCC [     3,    0,    0,    0,    0,    2,    1,    1 ]
    SOCC [     0,    0,    1,    1,    0,    0,    0,    0 ]

  Energy converged.

  @UHF Final Energy:  -147.63402658708560
```
```
ORBITAL ENERGY EIGENVALUES
[-20.44092476 -20.43982562  -1.61895614  -1.0904038   -0.71572756  -0.71572756  -0.60782104  -0.41366455  -0.41366455   0.68744582]
[-20.40935163 -20.40920055  -1.48824156  -0.89768363  -0.55493306  -0.44377123  -0.44377123   0.27095791   0.27095791   0.7772591 ]

ORBITAL COEFFICIENTS
[[ -7.02723743e-01  -7.03498943e-01   1.71186239e-01  -1.88587548e-01  -4.38136920e-16   1.19802374e-17  -6.20530364e-02  -8.32248930e-16   2.84706141e-18   8.93225221e-02]
 [ -2.07586114e-02  -1.33148014e-02  -5.81674920e-01   7.92133414e-01   1.29734380e-15   3.71511481e-17   3.16982472e-01   3.08736402e-15  -3.65588944e-16  -5.57900536e-01]
 [  1.59916666e-17  -2.52119849e-16   1.27385408e-15  -7.42459545e-16  -4.33303675e-01  -4.96366048e-01  -2.08215235e-15   7.60305796e-01  -1.06614927e-01   1.76041920e-16]
 [ -1.37080901e-16   3.00960019e-16  -1.03883038e-15   1.24001296e-16  -4.96366048e-01   4.33303675e-01  -1.09619443e-15   1.06614927e-01   7.60305796e-01  -5.66762209e-17]
 [  5.90696987e-03   3.16583297e-05   1.57503690e-01   1.31270800e-01  -2.26167345e-15   5.19105390e-16   6.18859802e-01   4.64771775e-16  -6.38429068e-16   9.47820940e-01]
 [  7.02723743e-01  -7.03498944e-01   1.71186239e-01   1.88587548e-01  -2.58412115e-16  -9.27419675e-17  -6.20530364e-02  -3.28901489e-16   9.31297183e-17  -8.93225221e-02]
 [  2.07586114e-02  -1.33148014e-02  -5.81674920e-01  -7.92133414e-01   3.84888710e-16   4.50797609e-16   3.16982472e-01   1.42655930e-15  -5.81164887e-16   5.57900536e-01]
 [  3.44654784e-18   1.38685853e-16  -3.47376755e-15   7.39837664e-16  -4.33303675e-01  -4.96366048e-01   7.58971316e-16  -7.60305796e-01   1.06614927e-01  -3.10964058e-16]
 [  5.02204089e-17  -2.41135394e-16  -1.66280656e-15   5.28441640e-16  -4.96366048e-01   4.33303675e-01  -2.13054056e-15  -1.06614927e-01  -7.60305796e-01  -1.82929277e-16]
 [  5.90696987e-03  -3.16583298e-05  -1.57503690e-01   1.31270800e-01   1.62284816e-15  -6.40098778e-16  -6.18859802e-01  -5.62572989e-16   7.22391062e-16   9.47820940e-01]]

[[  7.03199293e-01   7.03825068e-01   1.66962912e-01   1.81606050e-01  -6.94242840e-02   1.26199526e-16   2.37550839e-16   1.77535693e-16  -2.92964372e-16  -9.94707578e-02]
 [  1.85908180e-02   1.21320747e-02  -5.66324587e-01  -7.58247086e-01   3.43700682e-01  -6.01633883e-16  -1.20255279e-15  -9.26041850e-16   1.11962646e-15   6.03220533e-01]
 [ -2.30303447e-17   3.74738486e-17   1.00456372e-16   3.85334205e-16   2.57053817e-16  -6.56094837e-01   6.05879007e-02   3.82223195e-01   6.65835622e-01  -1.14402786e-15]
 [  1.36810067e-16  -6.76093540e-17  -6.51507432e-16   2.75583445e-16  -4.07776121e-15  -6.05879007e-02  -6.56094837e-01  -6.65835622e-01   3.82223195e-01  -9.99012108e-16]
 [ -5.34087079e-03  -1.80846940e-04   1.86114521e-01  -1.86377670e-01   6.10865141e-01   5.78752902e-16   1.76602239e-17  -2.60046091e-15   6.16760811e-16  -9.38544772e-01]
 [ -7.03199293e-01   7.03825068e-01   1.66962912e-01  -1.81606050e-01  -6.94242840e-02  -1.09924417e-16  -8.22394905e-17   1.92915095e-16  -2.94288502e-17   9.94707578e-02]
 [ -1.85908180e-02   1.21320747e-02  -5.66324587e-01   7.58247086e-01   3.43700682e-01   6.32390267e-16   4.46164071e-16  -9.94462356e-16   1.50639288e-16  -6.03220533e-01]
 [ -4.73773095e-18   6.29499495e-17   7.34517638e-18  -4.13141002e-17  -1.13153197e-16  -6.56094837e-01   6.05879007e-02  -3.82223195e-01  -6.65835622e-01   2.06338813e-16]
 [ -4.63093672e-17  -1.75133756e-17   6.65940253e-16   1.48690066e-16   2.63853965e-15  -6.05879007e-02  -6.56094837e-01   6.65835622e-01  -3.82223195e-01  -3.54232274e-16]
 [ -5.34087079e-03   1.80846940e-04  -1.86114521e-01  -1.86377670e-01  -6.10865141e-01   4.88140539e-16   1.51136604e-15   2.46856157e-15  -2.14894630e-15  -9.38544772e-01]]

TOTAL NUCLEAR REPULSION ENERGY: 28.0036243561 a.u.
TOTAL ENERGY: -147.634028138 a.u.

*****************************************************************************************************

TIME TAKEN: 7.575840318746571
```

### Density Functional Theory
##### SVWN He STO-3G
Completed DFT/SVWN calculation of He using the STO-3G basis set. Fairly close agreement to Psi4, differences are probably due to the different numerical integration techniques used.
```
  @RKS Final Energy:    -2.80959859524104

   => Energetics <=

    Nuclear Repulsion Energy =              0.0000000000000000
    One-Electron Energy =                  -3.8634969002750466
    Two-Electron Energy =                   2.1114258854701449
    DFT Exchange-Correlation Energy =      -1.0575275804361373
    Empirical Dispersion Energy =           0.0000000000000000
    PCM Polarization Energy =               0.0000000000000000
    EFP Energy =                            0.0000000000000000
    Total Energy =                         -2.8095985952410389
```
```
NUCLEAR REPULSION ENERGY:    0.0 a.u.
SCF ENERGY:                  -2.809831281891216 a.u.
CORRELATION ENERGY:          0.0 a.u.
TOTAL ENERGY:                -2.80983128189 a.u.

*************************************************************************************************

TIME TAKEN: 48.94459470800236s
```

##### SVWN H<sub>2</sub> STO-3G
```
  @RKS Final Energy:    -1.15582107524090

   => Energetics <=

    Nuclear Repulsion Energy =              0.7559674408428576
    One-Electron Energy =                  -2.5557060084773302
    Two-Electron Energy =                   1.3647790651283140
    DFT Exchange-Correlation Energy =      -0.7208615727347415
    Empirical Dispersion Energy =           0.0000000000000000
    PCM Polarization Energy =               0.0000000000000000
    EFP Energy =                            0.0000000000000000
    Total Energy =                         -1.1558210752409002
```
```
NUCLEAR REPULSION ENERGY:    0.755967561438 a.u.
SCF ENERGY:                  -1.9118948955020598 a.u.
CORRELATION ENERGY:          0.0 a.u.
TOTAL ENERGY:                -1.15592733406 a.u.

*************************************************************************************************

TIME TAKEN: 46.66810950479239s
```

### Møller–Plesset Second Order
##### CO STO-3G
MP2 calculation of CO with a bond-length of 2.14005 a<sub>0</sub> using the STO-3G basis set. Comparison with Psi4,
```
        Computing MP2 energy using SCF MOs (Canonical MP2)... 
        ============================================================================== 
        Nuclear Repulsion Energy (a.u.)    :    22.42938094724525
        SCF Energy (a.u.)                  :  -111.22496033969136
        REF Energy (a.u.)                  :  -111.22496033969136
        Alpha-Alpha Contribution (a.u.)    :    -0.01678813489592
        Alpha-Beta Contribution (a.u.)     :    -0.09597521248030
        Beta-Beta Contribution (a.u.)      :    -0.01678813489592
        Scaled_SS Correlation Energy (a.u.):    -0.01119208993061
        Scaled_OS Correlation Energy (a.u.):    -0.11517025497636
        SCS-MP2 Total Energy (a.u.)        :  -111.35132268459833
        SOS-MP2 Total Energy (a.u.)        :  -111.22496033969136
        SCSN-MP2 Total Energy (a.u.)       :  -111.28405457452502
        SCS-MP2-VDW Total Energy (a.u.)    :  -111.36459674656207
        SOS-PI-MP2 Total Energy (a.u.)     :  -111.35932563716378
        MP2 Correlation Energy (a.u.)      :    -0.12955148227214
        MP2 Total Energy (a.u.)            :  -111.35451182196350
        ============================================================================== 
```
```
NUCLEAR REPULSION ENERGY:    22.4293809472 a.u.
SCF ENERGY:                  -133.65434203636403 a.u.
CORRELATION ENERGY:          -0.129556206728 a.u.
TOTAL ENERGY:                -111.354517296 a.u.

*****************************************************************************************************

TIME TAKEN: 36.23057195376635s
```

### Coupled Cluster Singles and Doubles
##### CH<sub>4</sub> STO-3G
CCSD calculation of CH<sub>4</sub> in the STO-3G basis set. comparison with Psi4,
```
        SCF energy       (wfn)                =  -39.726835850063679
        Reference energy (file100)            =  -39.726835850063708

        Opposite-spin MP2 correlation energy  =   -0.053165399377273
        Same-spin MP2 correlation energy      =   -0.002972275669006
        MP2 correlation energy                =   -0.056137675046280
      * MP2 total energy                      =  -39.782973525109988

        Opposite-spin CCSD correlation energy =   -0.076647951103079
        Same-spin CCSD correlation energy     =   -0.001821943731474
        CCSD correlation energy               =   -0.078469894846414
      * CCSD total energy                     =  -39.805305744910122
```
```
NUCLEAR REPULSION ENERGY:    13.486321423 a.u.
SCF ENERGY:                  -53.21315722943018 a.u.
CORRELATION ENERGY:          -0.0784698838517 a.u.
TOTAL ENERGY:                -39.8053056903 a.u.

*************************************************************************************************

TIME TAKEN: 354.9196610947245s
```

## Citations
Spartan Student Wavefunction, Inc. Irvine, CA Except for molecular mechanics and semi-empirical models, the calculation methods used in Spartan Student have been documented in: Y. Shao, L.F. Molnar, Y. Jung, J. Kussmann, C. Ochsenfeld, S.T. Brown, A.T.B. Gilbert, L.V. Slipchenko, S.V. Levchenko, D.P. O’Neill, R.A. DiStasio Jr., R.C. Lochan, T. Wang, G.J.O. Beran, N.A. Besley, J.M. Herbert, C.Y. Lin, T. Van Voorhis, S.H. Chien, A. Sodt, R.P. Steele, V.A. Rassolov, P.E. Maslen, P.P. Korambath, R.D. Adamson, B. Austin, J. Baker, E.F.C. Byrd, H. Dachsel, R.J. Doerksen, A. Dreuw, B.D. Dunietz, A.D. Dutoi, T.R. Furlani, S.R. Gwaltney, A. Heyden, S. Hirata, C-P. Hsu, G. Kedziora, R.Z. Khalliulin, P. Klunzinger, A.M. Lee, M.S. Lee, W.Z. Liang, I. Lotan, N. Nair, B. Peters, E.I. Proynov, P.A. Pieniazek, Y.M. Rhee, J. Ritchie, E. Rosta, C.D. Sherrill, A.C. Simmonett, J.E. Subotnik, H.L. Woodcock III, W. Zhang, A.T. Bell, A.K. Chakraborty, D.M. Chipman, F.J. Keil, A.Warshel, W.J. Hehre, H.F. Schaefer, J. Kong, A.I. Krylov, P.M.W. Gill and M. Head-Gordon, Phys. Chem. Chem. Phys., 8, 3172 (2006).

“Psi4: An open-source ab initio electronic structure program,” J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein, F. Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke, M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl, W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill, and T. D. Crawford, WIREs Comput. Mol. Sci. 2, 556 (2012).