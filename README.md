##Introduction
A basic quantum chemical program written in python using the numpy and scipy libraries and is currently at the very beginnings. My time is tight as I try to balance this between my job and real-life so development may go slowly.

The basic HF part is complete, calculations can be carried out on all molecules and basis sets. It's a little slow at the moment and I need to work on making it faster. The biggest problem at the moment are the ERI, currently having a look at different method and trying to implement them. I have some plans to move the project to C++ but I think its better for me to understand computational chemistry first and then learn another language as python will speed up the process. 

I'm basing this work on Attlia Szabo and Neil S. Ostlunds "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory". It is a old but well written book, and is by far the best introduction to the computational chemistry I have read. The only other textbook that comes close to this is probably Robert G. Parr and Weitao Yangs "Density Functional Theory of Atoms and Molecules". I'm also looking at David B. Cooks "Handbook of Computational Chemistry" which appears to be quite in depth and has some programming parts.

The tests section is very poor at the moment because it isn't worth my time will start adding more later.

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
  
and from my program,
    
    TOTAL NUCLEAR REPULSION ENERGY: 1.36687066232 a.u.
    TOTAL ENERGY: -2.92922773418 a.u.
    
    *********************************************************************************************************
    
    Time Taken: 10.986138812057659s

##Calculation of C<sub>2</sub>H<sub>4</sub> in 3-21G
With some tweaking the code is fast enough to calculate C<sub>2</sub>H<sub>4</sub>.

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

Python Program,
    
    TOTAL NUCLEAR REPULSION ENERGY: 33.4424184132 a.u.
    TOTAL ENERGY: -77.6004608443 a.u.
    
    *********************************************************************************************************
    
    Time Taken: 72.70841306778921s

##UHF of O<sub>2</sub> in STO-3G
UHF calculation have been implemented, the below shows energy calculations of O<sub>2</sub> with a multiplicity of 3. When calculating in the Psi4,

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

and from Python,

    BEGIN SCF PROCEDURE
    SCF ENERGY: -174.075470966916 a.u.
    SCF ENERGY: -175.53766748081466 a.u.
    SCF ENERGY: -175.6225967630879 a.u.
    SCF ENERGY: -175.6352843627404 a.u.
    SCF ENERGY: -175.63727570713576 a.u.
    SCF ENERGY: -175.63759236091536 a.u.
    SCF ENERGY: -175.6376428896343 a.u.
    SCF ENERGY: -175.63765096020234 a.u.
    SCF ENERGY: -175.63765224957692 a.u.
    SCF ENERGY: -175.6376524555832 a.u.
    SCF ENERGY: -175.63765248849757 a.u.
    SCF ENERGY: -175.63765249375638 a.u.
    SCF ENERGY: -175.6376524945966 a.u.
    TIME TAKEN: 0.403388799599643s
    
    
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