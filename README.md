##Introduction
A basic quantum chemical program written in python using the numpy and scipy libraries and is currently at the very beginnings. My time is tight as I try to balance this between my job and real-life so development may go slowly.

The basic HF part is complete, calculations can be carried out on all molecules and basis sets. It's a little slow at the moment and I need to work on making it faster. The biggest problem at the moment are the ERI, currently having a look at different method and trying to implement them. I have some plans to move the project to C++ but I think its better for me to understand computational chemistry first and then learn another language as python will speed up the process. 

I'm basing this work on Attlia Szabo and Neil S. Ostlunds "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory". It is a old but well written book, and is by far the best introduction to the computational chemistry I have read. The only other textbook that comes close to this is probably Robert G. Parr and Weitao Yangs "Density Functional Theory of Atoms and Molecules". I'm also looking at David B. Cooks "Handbook of Computational Chemistry" which appears to be quite in depth and has some programming parts.

The tests section is very poor at the moment because it isn't worth my time will start adding more later.

Please see my GitHub pages http://chicheng45.github.io/Quantum_Chemistry/ for documentation and comparisons against other Qunatum Chemical Software Packages. 