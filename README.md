# Lid-Driven-Cavity
MATLAB code that solves the planar Navier-Stokes equations using the FMG-FAS algorithm of Smith (2023) for the Lid Driven Cavity problem.

The PDEs are in Streamfunction-Vorticity form, which decreases the dimension by one, but also makes boundary conditions difficult to implement (so some care and attention to detail is needed). These are then discretised via finite differences producing a large system of non-linear equations. As a full scale Newton method requires a huge amount of computational resources, iterative methods are needed. This code utilises a geometric multigrid method to solve large systems of non-linear equations in a reasonably quick manner by iterating on smaller/coarser grids and correcting the solutions on larger/finer grids. The full details are in Smith (2023).

To run the code in its current format simply download all the files into the same directory and then run _LDC(**Re**)_ where _**Re**_ is the Reynolds number. After the code runs, all the data will be saved to a file called _LDC_Re=**Re**.mat_. To view the results simply run _figures.m_. In order to test the consistency and accuracy of the code, some results from [Ghia _et al_ (1982)](https://doi.org/10.1016/0021-9991(82)90058-4) are attached in a spreadsheet.

<p align="center">
  <img height="400" src="https://github.com/bensmith95/Lid-Driven-Cavity/assets/29705711/b29c1253-2eb8-4716-b6eb-7f67e3310165">
</p>
