##Introduction
A basic quantum chemical program written in python using the numpy and scipy library and is currently at the very beginnings. My time is tight as I try to balance this between my job and real-life so development may go slowly.

Currently I'm looking to write HF in python and then to more accurate methods such as CI, CCSD and DFT. After this project is complete I may try out other languages such as Fortran, C, C++, Java, C# or even Haskell if I have the time.

I'm basing this work on Attlia Szabo and Neil S. Ostlunds "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory". It is a old but well written book, and is by far the best introduction to the computational chemistry I have read. The only other textbook that comes close to this is probably Robert G. Parr and Weitao Yangs "Density Functional Theory of Atoms and Molecules". After this I intend to check out David B. Cooks "Handbook of Computational Chemistry" which appears to be quite in depth.

I must apologies for the bad/ugly code, I'm going to try refactor the code sometime in future and begin developing properly when I get more comfortable with writing and unit testing in python.

##Calculation of HeH<sup>+</sup> in STO-3G
Completed my first calculation of HeH<sup>+</sup> with a bond-length of 1.4632 a<sub>0</sub> using the STO-3G basis set. Below is the output from Spartan Student Edition v5,

    SCF model:
     A restricted Hartree-Fock SCF calculation will be
     performed using Pulay DIIS + Geometric Direct Minimization
    
     SCF total energy:      -2.8418365 hartrees
    
      Reason for exit: Successful completion 
      Quantum Calculation CPU Time :          .39
      Quantum Calculation Wall Time:         2.85

and the output from my program,
  
    ORBITAL COEFFICIENTS
    [[-0.87660574  0.79774813]
     [-0.20247895 -1.16783645]]
    
    TOTAL ENERGY: -2.841836483 a.u.
    
    
    Time Taken: 0.10287352007652968s

##Calculation of HeH<sup>+</sup> in 3-21G

    SCF model:
     A restricted Hartree-Fock SCF calculation will be
     performed using Pulay DIIS + Geometric Direct Minimization
    
     SCF total energy:      -2.8873832 hartrees
    
      Reason for exit: Successful completion 
      Quantum Calculation CPU Time :          .58
      Quantum Calculation Wall Time:         1.28

and from my program,
    
    ORBITAL COEFFICIENTS
    [[-0.46642939  0.20194962 -0.02511278  1.1790321 ]
     [-0.52693222  0.77891794  0.0516036  -1.43875342]
     [-0.20335162 -0.33518808 -1.25188022 -0.06982107]
     [-0.0198215  -1.07250037  1.07697412  0.59541063]]
    
    SCF ENERGY: -4.254254337971538 a.u.
    TOTAL NUCLEAR REPULSION ENERGY: 1.36686714051 a.u.
    TOTAL ENERGY: -2.88738719746 a.u.

##Calculation of HeH<sup>+</sup> in 6-311+G**
Able to run calculations for all types of orbitals now, below shows the HeH+ molecule using the 6-311+G basis set which contains polarization functions for H and He. This calculation in particular takes quite a lot longer then Spartan and has a lot of inefficiencies that need sorting out. 
        
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
    
      Reason for exit: Successful completion 
      Quantum Calculation CPU Time :          .36
      Quantum Calculation Wall Time:          .42
    
    SPARTAN STUDENT Properties Program:  (PC/x86)                  Release  5.0.0  
      Use of molecular symmetry enabled
    
                         Cartesian Coordinates (Angstroms)
           Atom            X             Y             Z     
        ---------    ------------- ------------- -------------
    
      1 He He1          0.0000000     0.0000000     0.2580967
      2 H  H2           0.0000000     0.0000000    -0.5161934
    
      Point Group = CIV Order =  1 Nsymop =  5
    
      Closed-Shell Molecular Orbital Coefficients
      MO:                   1          2          3          4          5    
      Eigenvalues:      -1.63348   -0.26340    0.00717    0.51943    0.73352
              (ev)     -44.44922   -7.16735    0.19506   14.13428   19.96021
    
                           A1         A1         A1         A1         E1x    
       1 He1   S         0.25858   -0.11142   -0.00273   -0.15042    0.00000
       2 He1   S'        0.48928   -0.21808    0.02545   -1.02708    0.00000
       3 He1   S"        0.23234   -0.77157    0.37426    2.19847    0.00000
       4 He1   PX        0.00000    0.00000    0.00000    0.00000    0.61780
       5 He1   PY        0.00000    0.00000    0.00000    0.00000    0.00000
       6 He1   PZ       -0.03666   -0.01140   -0.02567   -0.18596    0.00000
       7 H2    S         0.12788    0.17096   -0.13705    0.08557    0.00000
       8 H2    S'        0.11737    0.59194   -1.39013   -0.17403    0.00000
       9 H2    S"       -0.01494    0.81277    1.34904   -1.00788    0.00000
      10 H2    PX        0.00000    0.00000    0.00000    0.00000    0.55682
      11 H2    PY        0.00000    0.00000    0.00000    0.00000    0.00000
      12 H2    PZ        0.03386   -0.00708   -0.08120   -0.27526    0.00000
      MO:                   6          7          8          9         10    
      Eigenvalues:       0.73352    1.24984    1.49540    1.49540    1.79599
              (ev)      19.96021   34.00979   40.69199   40.69199   48.87150
    
                           E1y        A1         E1x        E1y        A1     
       1 He1   S         0.00000    0.03722    0.00000    0.00000   -0.00357
       2 He1   S'        0.00000    0.14439    0.00000    0.00000   -1.11642
       3 He1   S"        0.00000   -0.53891    0.00000    0.00000   -0.54953
       4 He1   PX        0.00000    0.00000   -0.93246    0.00000    0.00000
       5 He1   PY        0.61780    0.00000    0.00000   -0.93246    0.00000
       6 He1   PZ        0.00000    0.97627    0.00000    0.00000    0.73179
       7 H2    S         0.00000   -0.37178    0.00000    0.00000   -0.99342
       8 H2    S'        0.00000    1.50195    0.00000    0.00000    2.61899
       9 H2    S"        0.00000   -0.52957    0.00000    0.00000   -0.66424
      10 H2    PX        0.00000    0.00000    0.97011    0.00000    0.00000
      11 H2    PY        0.55682    0.00000    0.00000    0.97011    0.00000
      12 H2    PZ        0.00000   -0.24309    0.00000    0.00000    1.41949
      MO:                  11         12    
      Eigenvalues:       2.63499    6.32407
              (ev)      71.70167  172.08679
    
                           A1         A1     
       1 He1   S        -0.15831    1.50919
       2 He1   S'       -0.83288   -2.39836
       3 He1   S"       -1.23816    0.54728
       4 He1   PX        0.00000    0.00000
       5 He1   PY        0.00000    0.00000
       6 He1   PZ        1.76790    0.65768
       7 H2    S         1.21322    0.33605
       8 H2    S'        0.96289    0.75430
       9 H2    S"        0.28477   -0.34835
      10 H2    PX        0.00000    0.00000
      11 H2    PY        0.00000    0.00000
      12 H2    PZ        1.62046    0.74519
       
       
    
      Reason for exit: Successful completion 
      Properties CPU Time :          .23
      Properties Wall Time:          .21
  
and from my program,

    *********************************************************************************************************
    
    BEGIN SCF PROCEDURE
    SCF ENERGY: -4.197619917017029 a.u.
    SCF ENERGY: -4.295101268737092 a.u.
    SCF ENERGY: -4.2960790541862846 a.u.
    SCF ENERGY: -4.296097875302424 a.u.
    SCF ENERGY: -4.29609837993943 a.u.
    SCF ENERGY: -4.296098395922256 a.u.
    SCF ENERGY: -4.296098396493559 a.u.
    
    ORBITAL ENERGY EIGENVALUES
    [[-1.63347727  0.          0.          0.          0.          0.          0.          0.          0.          0.          0.          0.        ]
     [ 0.         -0.26339532  0.          0.          0.          0.          0.          0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.00716821  0.          0.          0.          0.          0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.19356385  0.          0.          0.          0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.51942552  0.          0.          0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.99158389  0.          0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.          1.24983604  0.          0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.          0.          1.27372407  0.          0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.          0.          0.          1.79599137  0.          0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.          0.          0.          0.          1.99898165  0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.          0.          0.          0.          0.          2.63498222  0.        ]
     [ 0.          0.          0.          0.          0.          0.          0.          0.          0.          0.          0.          6.32406869]]
    
    ORBITAL COEFFICIENTS
    [[  2.58584881e-01   1.11419466e-01   2.73005263e-03  -7.48702927e-19  -1.50423550e-01  -2.17948523e-16  -3.72200990e-02  -2.43457458e-16  -3.57433126e-03   4.31291086e-17  -1.58309597e-01   1.50918881e+00]
     [  4.89277720e-01   2.18077259e-01  -2.54503646e-02  -3.89395314e-17  -1.02707500e+00  -1.71143424e-15  -1.44392943e-01  -1.33010178e-15  -1.11641553e+00  -1.07023484e-14  -8.32880510e-01  -2.39835532e+00]
     [  2.32338183e-01   7.71565750e-01  -3.74263305e-01   8.12637263e-17   2.19846827e+00   1.30572863e-15   5.38904289e-01   4.02502750e-15  -5.49526145e-01  -5.29239641e-15  -1.23816513e+00   5.47280694e-01]
     [ -1.25675432e-16   2.56879167e-17  -1.04247031e-16   8.03584886e-02  -9.91392876e-17  -5.13306318e-01   5.08038445e-15  -7.69007413e-01  -5.63271700e-15   8.22484898e-01   2.64406368e-16  -4.66076768e-16]
     [  6.53480280e-17  -1.62541951e-16  -3.27766428e-17   3.00963956e-01   1.68977446e-16  -5.92396533e-01  -1.27676743e-15   2.53213997e-01   5.17763633e-15  -6.77431289e-01   1.40030086e-15   1.78471452e-16]
     [ -3.66592789e-02   1.14004491e-02   2.56705689e-02   2.97908164e-17  -1.85959507e-01   9.46754943e-16  -9.76269455e-01  -7.36822876e-15   7.31784078e-01   7.11763402e-15   1.76790563e+00   6.57679110e-01]
     [  1.27877984e-01  -1.70963252e-01   1.37050368e-01  -6.19877633e-17   8.55676316e-02   5.03408197e-16   3.71779872e-01   2.10088016e-15  -9.93418756e-01  -8.23540608e-15   1.21321808e+00   3.36052584e-01]
     [  1.17371441e-01  -5.91943042e-01   1.39013492e+00   7.69131846e-17  -1.74032110e-01   1.63879483e-15  -1.50194617e+00  -1.04975903e-14   2.61898375e+00   2.31358612e-14   9.62907113e-01   7.54303786e-01]
     [ -1.49407535e-02  -8.12768101e-01  -1.34903802e+00  -8.21709122e-17  -1.00787879e+00  -1.70029375e-15   5.29569491e-01   3.87901327e-15  -6.64240646e-01  -5.17130086e-15   2.84762775e-01  -3.48352445e-01]
     [  2.17834896e-16   1.95532778e-18   1.01103199e-16   7.45884503e-02   7.42025366e-17   4.58927764e-01   5.41556012e-15  -7.25045215e-01   7.04395046e-15  -8.92146837e-01  -8.39941456e-16   5.15019913e-16]
     [ -7.87827883e-17   1.66899763e-16   3.15134977e-17   2.82351090e-01  -2.07608673e-16   6.14132132e-01  -1.37417874e-15   1.43240922e-01  -5.75952414e-15   6.97957777e-01  -1.10206255e-15  -1.73449863e-16]
     [  3.38599636e-02   7.08226655e-03   8.12018394e-02  -1.33539218e-17  -2.75256643e-01   1.28549614e-15   2.43088515e-01   2.02747055e-15   1.41948897e+00   1.33137247e-14   1.62046309e+00   7.45186953e-01]]
    
    TOTAL NUCLEAR REPULSION ENERGY: 1.36687066232 a.u.
    TOTAL ENERGY: -2.92922773418 a.u.
    
    *********************************************************************************************************
    
    Time Taken: 10.986138812057659s