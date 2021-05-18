## Introduction
A basic quantum chemistry program written in Python3 using the numpy and scipy libraries.

Currently this program fully supports RHF, UHF, CIS, TDHF, DFT, CCSD and CCSD(T). My next plans are to implement more DFT functionals.

I used Attlia Szabo and Neil S. Ostlunds _Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory_ and David B. Cooks _Handbook of Computational Quantum Chemistry_ as my main references for the theories and methods behind the electronic structure calculations. The developers resources at http://www.psicode.org/developers.php were also invaluable to the success of project and had a number of excellent tutorials and programming examples.

## Instructions
Install the Gaussium package using
```text
python setup.py install
```
you can then run Gaussium using the start function for example
```python
from gaussium import start


def run():
    start('He.mol', '3-21G.gbs', ('DFT', 'S', 'VWN5'), 4)  # -2.80601675458 a.u.


if __name__ == "__main__":
    run()
``` 
assuming that the molecule and basis set files are in the working directory. This calculation will start a DFT calculation with the SVWN5 functional with the 3-21G basis set for a Helium atom and will use four processors, see the examples directory for more examples.

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
