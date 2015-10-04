##Introduction
A basic quantum chemical program written in python using the numpy and scipy library and is currently at the very beginnings. My time is tight as I try to balance this between my job and real-life so development may go slowly.

Currently I'm looking to write HF in python and then to more accurate methods such as CI, CCSD and DFT. After this project is complete I may try out other languages such as Fortran, C, C++, Java, C# or even Haskell if I have the time.

I'm basing this work on Attlia Szabo and Neil S. Ostlunds "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory". It is a old but well written book, and is by far the best introduction to the computational chemistry I have read. The only other textbook that comes close to this is probably Robert G. Parr and Weitao Yangs "Density Functional Theory of Atoms and Molecules". After this I intend to check out David B. Cooks "Handbook of Computational Chemistry" which appears to be quite in depth.

I must apologies for the bad/ugly code, I'm going to try refactor the code sometime in future and begin developing properly when I get more comfortable with writing and unit testing in python.

##Calculation of HeH<sup>+</sup> in STO-3G
Completed my first calculation of HeH<sup>+</sup> with a bond-length of 1.4632 a<sub>0</sub> using the STO-3G basis set. Below is the output from Spartan Student Edition v5,

    SPARTAN STUDENT Quantum Mechanics Program:  (PC/x86)   Release  5.0.0v4
    
    Job type: Single point.
    Method: RHF
    Basis set: STO-3G
    Number of shells: 2
    Number of basis functions: 2
    Charge :     +1 
    Multiplicity: 1
    
    SCF model:
     A restricted Hartree-Fock SCF calculation will be
     performed using Pulay DIIS + Geometric Direct Minimization
    
     SCF total energy:      -2.8418365 hartrees
    
      Reason for exit: Successful completion 
      Quantum Calculation CPU Time :          .39
      Quantum Calculation Wall Time:         2.85
    
    SPARTAN STUDENT Properties Program:  (PC/x86)                  Release  5.0.0  
       
       
    
      Reason for exit: Successful completion 
      Properties CPU Time :          .28
      Properties Wall Time:          .37

and the output from my program,

    ----------ITERATION: 6----------
    
    DENSITY MATRIX
    [[ 1.53690875  0.35496777]
     [ 0.35496777  0.08198412]]
    
    FOCK MATRIX
    [[-1.59018101 -1.06102034]
     [-1.06102034 -0.83401647]]
    
    ORBITAL ENERGY EIGENVALUES
    [[-1.63279763  0.        ]
     [ 0.         -0.17248459]]
    
    ORBITAL COEFFICIENTS
    [[-0.87660574  0.79774813]
     [-0.20247895 -1.16783645]]
    
    TOTAL ENERGY: -2.841836483 a.u.
    
    
    Time Taken: 0.10287352007652968s

##Calculation of HeH<sup>+</sup> in 3-21G

    SPARTAN STUDENT Quantum Mechanics Program:  (PC/x86)   Release  5.0.0v4
    
    Job type: Reading previous wavefunction
    
    Job type: Single point.
    Method: RHF
    Basis set: 3-21G(*)
    Number of shells: 4
    Number of basis functions: 4
    Charge :     +1 
    Multiplicity: 1
    
    SCF model:
     A restricted Hartree-Fock SCF calculation will be
     performed using Pulay DIIS + Geometric Direct Minimization
    
     SCF total energy:      -2.8873832 hartrees
    
      Reason for exit: Successful completion 
      Quantum Calculation CPU Time :          .58
      Quantum Calculation Wall Time:         1.28
    
    SPARTAN STUDENT Properties Program:  (PC/x86)                  Release  5.0.0  
      Use of molecular symmetry enabled
    
                         Cartesian Coordinates (Angstroms)
           Atom            X             Y             Z     
        ---------    ------------- ------------- -------------
    
      1 He He1          0.0000000     0.0000000     0.2580001
      2 H  H2           0.0000000     0.0000000    -0.5160002
    
      Point Group = CIV Order =  1 Nsymop =  5
    
      Closed-Shell Molecular Orbital Coefficients
      MO:                   1          2          3          4    
      Eigenvalues:      -1.62778   -0.25548    0.56417    1.70220
              (ev)     -44.29414   -6.95208   15.35197   46.31923
    
                           A1         A1         A1         A1     
       1 He1   S         0.46645   -0.20197    0.02520    1.17905
       2 He1   S'        0.52685   -0.77927   -0.05218   -1.43896
       3 H2    S         0.20345    0.33512    1.25189   -0.07012
       4 H2    S'        0.01976    1.07284   -1.07663    0.59585
       
       
    
      Reason for exit: Successful completion 
      Properties CPU Time :          .23
      Properties Wall Time:          .22
      
and from my program
    
        ----------ITERATION: 9----------
    
    DENSITY MATRIX
    [[ 0.43514231  0.49149589  0.18979481  0.01843242]
     [ 0.49149589  0.55514761  0.21437439  0.02081953]
     [ 0.18979481  0.21437439  0.08278227  0.00803962]
     [ 0.01843242  0.02081953  0.00803962  0.00078079]]
    
    FOCK MATRIX
    [[-0.76840144 -1.52102481 -0.84558193 -0.91315206]
     [-1.52102481 -1.13668518 -0.82954991 -0.85510639]
     [-0.84558193 -0.82954991 -0.43701373 -0.7688132 ]
     [-0.91315206 -0.85510639 -0.7688132  -0.71802821]]
    
    ORBITAL ENERGY EIGENVALUES
    [[-1.62777906  0.          0.          0.        ]
     [ 0.          1.70220407  0.          0.        ]
     [ 0.          0.         -0.25548372  0.        ]
     [ 0.          0.          0.          0.56417493]]
    
    ORBITAL COEFFICIENTS
    [[-0.46644524  1.17904612  0.20196503 -0.02519964]
     [-0.52685263 -1.438957    0.77926896  0.0521786 ]
     [-0.20344813 -0.07012054 -0.33512332 -1.25189158]
     [-0.01975852  0.5958455  -1.07284528  1.0766338 ]]
    
    ELECTRON ENERGY: -4.254767212271082 a.u.
    TOTAL ENERGY: -2.88738421923 a.u.