## Introduction
A basic quantum chemical program written in Python3 using the numpy and scipy libraries.

Currently this program fully supports RHF, UHF, CIS, TDHF, DFT, CCSD and CCSD(T). My next plans are to implement DFT into the program.

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
