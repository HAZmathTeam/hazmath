# HAZMATH: A Simple Finite Element, Graph, and Solver Library

---

**Authors:** [Xiaozhe \***H**\*u (Tufts)](http://math.tufts.edu/faculty/xhu/), [James \***A**\*dler (Tufts)](http://math.tufts.edu/faculty/jadler), [Ludmil \***Z**\*ikatanov (Penn State)](http://personal.psu.edu/ltz1/)

## Contributions:

- ***HAZNICS (HAZMATH+FEniCS) and Python interface:*** Ana Budisa (Simula, Norway), Miroslav Kuchta (Simula, Norway), Kent-Andre Mardal (Simula, Univ Oslo, Norway).
- ***Rational Approximation of Functions***: Clemens Hofreither (RICAM, Austrian Academy of Sciences)
- ***Grid refinement***: Yuwen Li (Penn State)
- ***Geometric MultiGrid:*** Johannes Kraus (Universitat Duisburg-Essen, Germany), Peter Ohm (Tufts), Yunrong Zhu (Idaho State).

---

**Overview:** The HAZMATH Finite Element (FE) and Solver library is built from C (C99 compatible) source files and provides software components that can be used to simulate physical/social phenomena described by Partial Differential Equations (PDEs), systems of PDEs, or graphs. We have included several standard discretizations for PDEs, often used in various applications. Our examples range from simple scalar elliptic PDEs to systems of PDEs, include time-dependent and nonlinear problems, and contain methods for the solution of linear systems. 
 
Our aim is to provide a basic tool, which can be used to tackle specific problems on demand. We do not intend to create a universal package that contains all discretizations available. HAZMATH contains the software to solve numerical models based on finite-element/volume/difference discretizations, and our team can assist in tweaking and adjusting these basic tools to meet the demands of the specific application. Thus, we aim to complement the user's expertise and help build application-specific packages.

---

**Components:** As of early 2017, the main components of the HAZMATH FE/Graph library are:

1. **Basic FE:**  Low order (up to order 2) continuous FE and mixed FE discretizations for Darcy's flow and scalar elliptic equations; Stable discretizations for Stokes and linear elasticity; Discretization for full Maxwell's equations using Nedelec and Raviart-Thomas elements.

2. **Solvers:** Unsmoothed Aggregation Algebraic Multigrid Methods (UA-AMG); Preconditioned Krylov subspace methods.

3. **HAZNICS:** Collaboration with Simula; Interface with Python and use of HAZMATH solvers in applications modeled by [FEniCS](https://fenicsproject.org/) finite element libraries.

3. **Interfaces with External Libraries:** Some routines require a direct linear solver and these depend on the [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) library (by Tim Davis, Texas A & M), which the user needs to install. We also provide an interface for using [Multigraph 2.1](http://ccom.ucsd.edu/~reb/software.html) (by Randolph E Bank, UCSD) as a solver.  


---

**Obtaining the Software:** A tarball with installation instructions and the source files needed to build the HAZmath library and the examples can be requested by emailing the developers at  _hazmath [at] googlegroups [dot] com_. 


**Licensing:** This software is a free software distributed under the GNU's Lesser
General Public License or LGPL, Version 3 or any later versions. This
software is distributed WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the [GNU's Lesser General Public License](https://www.gnu.org/licenses/#LGPL) for more details.

**System Requirements:** The library should build on any standard Linux OS or MAC OS X, using cmake (>=3.12).

---

**TODO:** More detailed documentation; Standartize mesh formats; Minimum Stretching Tree Multigrid; Signal Processing on Graphs; 

---

(c) 2009- by X. Hu, J. Adler, L. Zikatanov 

---
