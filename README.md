# DGVlib
DGVlib is a collections of subroutines to enable nodal-DG discretization of the Boltzmann equation

Contributors: Alex Aleskeenko, Jeffrey Limbacher, Truong Nguyen, Craig Euler, Amy Grandilli,

E-mail: alexander.alekseenko@csun.edu

This is a copy of a Fortran code that implements a collection of tools for solution of the Boltzmann equation. The key 
approach is the using nodal-Discontinuous Galerkin discretization in velocity variable. This library has been used in a number of 
0D, 1d, and 2D solvers. The code includes several models for evaluating the Boltzmann collision operator: the full Boltzmann 
collision operator for binary collisions (hard--spheres are currently supported), linearized collision operator, BGk, ES-BGK, and 
model collision operator with enforced relaxation rates for moments. There are libraries used in modeling distribhution function as 
sums of maxwelliand and gaussians. Some subroutines were implemented as MPI parallel, some openMP. The goal is to have the most complehensive collection. 
Depending on how you wan tto use it, you may will not need all modules. 

Methods for discretization of the Boltzmann collision operator that are used in this code were reported in following papers:

2020 A. Alekseenko, A. Grandilli, A. Wood, An ultra-sparse approximation of kinetic solutions to spatially homogeneous flows of non-continuum gas, Results in Applied Mathematics, Vol. 5, 2020, 100085, ISSN 2590-0374, https://doi.org/10.1016/j.rinam.2019.100085.
2019 A. Alekseenko, J. Limbacher, Evaluating high order discontinuous Galerkin discretization of the Boltzmann collision integral in $O(N^2)$ operations using the discrete Fourier transform.\ \textit{Kinetic and Related Models} Vol.~12, n.~4,703-- 726, DOI: 10.3934/krm.2019027
2018 A.\ Alekseenko, T.\ Nguyen, and A.\ Wood, A deterministic-stochastic method for computing the Boltzmann collision Integral in $O(MN)$ operations. \textit{Kinetic and Related Models}, Vol., 11, no. 1937-5093
2015 A.\ Alekseenko and C.\ Euler, A Bhatnagar-Gross-Krook kinetic model with velocity-dependent collision frequency and corrected relaxation of moments.\textit{Continuum Mechanics and Thermodynamics} DOI: 10.1007/s00161-014-0407-0
2014 A.\ Alekseenko and E.\ Josyula, Deterministic solution of the spatially homogeneous Boltzmann equation using discontinuous Galerkin discretizations in the velocity space, \textit{Journal of Computational Physics}, Vol.~272, n.~1, (2014) 170--188.
2012 A.\ Alekseenko and E.\ Josyula, Deterministic solution of the Boltzmann equation using a discontinuous Galerkin velocity discretization.\ \textit{Proceedings of the $28^{\mathrm{th}}$ International Symposium on Rarefied Gas Dynamics, Spain 2012}, AIP Conference Proceedings, 2012, 8 pp.
