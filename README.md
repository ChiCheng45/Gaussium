##Introduction
A basic quantum chemical program written in python using the numpy and scipy library and is currently at the very beginnings. My time is tight as I try to balance this between my job and real-life so development may go slowly.

Currently I'm looking to write HF in python and then to more accurate methods such as CI, CCSD and DFT. After this project is complete I may try out other languages such as Fortran, C, C++, Java, C# or even Haskell if I have the time.

I'm basing this work on Attlia Szabo and Neil S. Ostlunds "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory". It is a old but well written book, and is by far the best introduction to the computational chemistry I have read. The only other textbook that comes close to this is probably Robert G. Parr and Weitao Yangs "Density Functional Theory of Atoms and Molecules". After this I intend to check out David B. Cooks "Handbook of Computational Chemistry" which appears to be quite in depth.

I must apologies for the bad/ugly code, I'm going to try refactor the code sometime in future and begin developing properly when I get more comfortable with writing and unit testing in python.

##Calculation of HeH+
Completed my first calculation of HeH+ with a bond-length of 1.4632 a_0 using the STO-3G basis set. Below is the output from Spartan Student Edition v5,

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
    
    Orbital Coefficients
    [[-0.87660574  0.79774813]
     [-0.20247895 -1.16783645]]
    
    TOTAL ENERGY: -2.841836483 a.u.
    
    
    Time Taken: 0.10287352007652968s