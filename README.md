# HAZMATH: A Simple Finite Element, Graph, and Solver Library

**Authors:** [Xiaozhe \***H**\*u (Tufts)](https://xiaozhehu.math.tufts.edu/), [James \***A**\*dler (Tufts)](https://jadler.math.tufts.edu/), [Ludmil \***Z**\*ikatanov (Penn State)](http://personal.psu.edu/ltz1/)

**Contributors:**

- ***HAZNICS (HAZMATH+FEniCS) and Python interface:*** Ana Budisa (Simula, Norway), Miroslav Kuchta (Simula, Norway), Kent-Andre Mardal (Simula, Univ Oslo, Norway).
- ***Rational Approximation of Functions***: Clemens Hofreither (RICAM, Austrian Academy of Sciences)
- ***Grid refinement and adaptive FE***: Yuwen Li (Penn State)
- ***Geometric MultiGrid:*** Johannes Kraus (Universitat Duisburg-Essen, Germany), Peter Ohm (Sandia), Yunrong Zhu (Idaho State).

---

**Overview:** The HAZMATH Finite Element (FE) and Solver library is built from C (C99 compatible) source files and provides software components that can be used to simulate physical/social phenomena described by Partial Differential Equations (PDEs), systems of PDEs, or graphs. We have included several standard discretizations for PDEs, often used in various applications. Our examples range from simple scalar elliptic PDEs to systems of PDEs, include time-dependent and nonlinear problems, and contain methods for the solution of linear systems. 
 
Our aim is to provide a basic tool, which can be used to tackle specific problems on demand. We do not intend to create a universal package that contains all discretizations available. HAZMATH contains the software to solve numerical models based on finite-element/volume/difference discretizations, and our team can assist in tweaking and adjusting these basic tools to meet the demands of the specific application. Thus, we aim to complement the user's expertise and help build application-specific packages.

---

**Components:** As of early 2017, the main components of the HAZMATH FE/Graph library are:

1. **Basic FE:**  Low order (up to order 2) continuous FE and mixed FE discretizations for Darcy's flow and scalar elliptic equations; Stable discretizations for Stokes and linear elasticity; Discretization for full Maxwell's equations using Nedelec and Raviart-Thomas elements.

2. **Solvers:** Unsmoothed Aggregation Algebraic Multigrid Methods (UA-AMG); Preconditioned Krylov subspace methods.

3. **HAZNICS:** Collaboration with Simula; Interface with Python and use of HAZMATH solvers in applications modeled by [FEniCS](https://fenicsproject.org/) finite element libraries.

4. **Interfaces with External Libraries:** Some routines require a direct linear solver and these depend on the [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) library (by Tim Davis, Texas A & M), which the user needs to install. We also provide an interface for using [Multigraph 2.1](http://ccom.ucsd.edu/~reb/software.html) (by Randolph E Bank, UCSD) as a solver.  

---

**System Requirements:** The library should build on any standard Linux OS or MAC OS X, using cmake (>=3.12).  After downloading type ```make help``` for a quick guide to building:

```
$$ make help
=============================================================
                         HAZMATH
-------------------------------------------------------------
Quick start: 1. Edit "haz_config/hazmath.mk" as needed.
             2. Type "make config" ;
                and then  "make install"
=============================================================
 make config           # set environment to build
 make config CC=gcc    # set a specific C compiler
 make config debug=yes # compiler's debug options ON
 make headers          # function defs in include/hazmath.h
 make docs	       # install documentations for HAZMATH
      		       # (!need to set doxygen=yes!)
 make                  # compile/link libhazmath
 make install          # install libhazmath and examples
 make clean            # remove *.o and *~ files
 make distclean        # delete the build directory
=============================================================
*More user defined options can go in "haz_config/hazmath.mk"*
```

---

**Licensing:** This software is a free software distributed under the MIT License. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. See further details about the [MIT License](https://opensource.org/licenses/MIT) at the [Open Source Initiative Licenses](https://opensource.org/licenses/) page.


(c) 2009- by X. Hu, J. Adler, L. Zikatanov 

---
