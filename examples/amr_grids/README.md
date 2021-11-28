
# INPUT FILES for the simple grid generator in HAZMATH

---

HAZMATH: A Simple Finite Element, Graph, and Solver Library

Copyright (c) 2009- HAZMath: Xiaozhe Hu, James H Adler, Ludmil T Zikatanov  

---

A short description of the input required to use the simple grid generator included with HAZMATH. Please refer to the examples for more details. 

## Rules ##

1. All empty lines are ignored, consecutive spaces and tabs are replaced by one space everything after a percent sign "%"  all the way to the end of the line is ignored.
2. Every input variable should be set as follows: **variable{value}**
3. If the variable is an array, the rule is the same (but **SPACE** or **TAB** separated values): **variable{value1 value2 ...}** 
4. If the input file is too long it is truncated to fixed number of bits (but there are enough bits to do very complicated meshes). 

Below we describe each group of parameters using as example an input file for a square domain in 2D. Many other examples are found in the input files in examples/amr_grids directory.

## Common parameters ##

**title{L-shaped domain}** (title)  
**dimension{2}**  (spatial dimension in which the grid has to be generated);  
**print_level{1}** (how much information to print on stdout (=0 less)||(>0 more));  
**refinement_type{0}** (non-negative means a criss cross grid; if negative, diagonals go on one side);  

### Filenames for output:##
As given below, also names can have spaces, etc.  
**file_grid{grid vtu/2d_L.grid}** (file where to write the resulting grid);  
**file_vtu{grid vtu/2dL.vtu}** (file where to write the vtu data for paraview).  

## AMR parameters ##

**amr_marking_type{1}** (how to mark for refinement).   
If {0}, in the example we presented, every simplex is refined; if non-zero,then only marked simplices are refined). In the example included with hazmath amr_marking{1} marks for refinement every simplex whose closure contains the origin. 

   
**num_refinements{16}** (how many refinements to make)  
**err_stop_refinement{-1.e-10}** (user defined parameter: stopping of the refinement using error estimator. Not used in the examples for now).

## Coordinate systems ##
**num_coordsystems{3}** (number of coordinate systems)  

 * Each coordinate system is represented by 2 integers and d real numbers in Rd. Below we have an example with one cartesian and two polar systems. The first polar is centered at (-1,-1) and the second at (0,0). In short:  
**coord_system**=(label  type  origin)

**data_coordsystems**  
**{0 0 0. 0.**  
**1  1 -1. -1.**  
**2  1 0. 0.}**
  
## Vertices ##

  * Every vertex is described by its number, its coord_system and its coordinates in this system; *vertex*=(label  coord_system  coords) as given below:  

**num_vertices{8}** (number of vertices describing the domain)  
**data_vertices**
**{0 0 -1.   -1.**  
** 1 0 -1.   0.**  
** 2 0  0.   -1.**  
** 3 0  0.   0.**  
** 4 0  1.   -1.**  	
** 5 0  1.    0.**  
** 6 0 -1.   1.**  
** 7 0  0.   1.}**


## Edges ##
  *  Edges and number of divisions for every edge. Not all edges need to be present in this list below. Parallel (after mapping to the unit cube) edges are divided by the same number of divisions to generate the mesh. Here if we have different number of divisions on parallel (after mapping to the unit cube) edges then the max number of divisions is taken. Also, if an edge is listed but does not exist, then it is ignored. If an edge is not present its divisions are set to 1. In short: edge=(v w num_divisions)  

**num_edges{6}**  (number of edges with prescribed divisions)  
**data_edges**  
**{0 1 3** %(v1 v2 num_divisions(edge))  
** 1 3 4**  
** 0 2 5**  
** 2 3 3**  
** 3 7 4**  
** 4 5 5**  
** 6 7 3}**

## Macroelements##

Macroelements are domains which are isomorphic to the unit cube via a d-linear transformation with non-zero Jaconian. Each macroelement has 2^{d} vertices.  **macroelement**=(v[0] ... v[2^{d}-1]  code(element))  
**num_macroelements{3}**  
**data_macroelements{2 4 3 5 1**  
** 0 1 2 3 2**  
** 6 7 1 3 -8}**

## Macrofaces ##

Similarly to macroelements, macrofaces are some (or all) of the the (d-1) dimensional faces of the macroelements. The boundary faces are detected automatically. **macroface**=(v[0] ... v[2^{d-1}-1]  code(face)). Faces that are listed here, but actually do not exist are ignored. Faces of macroelements that are not listed here are given  code 0 in the interior or code 1 (Dirichlet) on the boundary.

**num_macrofaces{10}** (number of macrofaces)  
  
**data_macrofaces**  
**{0 2  1** % face with code 1  
** 1 4  0** % face with code 0  
** 3 2  4**  
** 2 3  0**  
** 0 1  32**  
** 1 4  24**  
** 4 5  17**  
** 3 5  1**  
** 6 7  2**  
** 7 3  3}**


---

** Further details (skip perhaps...):**

NUMBERING OF VERTICES DESCRIBING A MACROELEMENT: The vertices should
be ordered so that the macroelement (polyhedron) is mapped to the unit
cube without singularities. This can be done by labeling all vertices
so that each face can be mapped to a face of the unit cube. In
dimension d, let us say that the labels of the macroelement vertices
are: (l[0]...l[(2^d-1)]). Such macroelement is mapped to the cube with
coordinates the binary representation of [0...(2^d-1)] i.e. in
coordinates, (x[0]...x[2^d-1])-->[(0,..,0)...(1,...,1)]. We call the
numbering of the macroelement vertices consistent if each
k-dimansional face of the macroelement is mapped to a k-dimensional
face in the unit cube for k<=d. Note that a k-dimensional face in the
unit cube consists of the vertices for which (d-k) bits in the binary
representation are fixed, i.e. (d-k) coordinates in the unit cube are
fixed.

For example, in 3D if a macroelement vertices are labeled as
(l[0]...l[7]), or, in binary (l[000]...l[111]). Then, fixing the
first bit to be 0 gives the face (l[0],l[1],l[2],l[3]), while
fixing the second bit to be 1 gives (l[2],l[3],l[6],l[7]), and so
on. The latter face is mapped to the back face of the unit cube,
i.e.  (x[l[2]],x[l[3]],x[l[6]],x[l[7]]) maps to
[2,3,6,7]=[010,011,110,111]=[(0,1,0),(0,1,1),(1,1,0),(1,1,1)] and
such mapping should also be without singlarity and this gives a
restriction on the numbering of the vertices because their
coordinates x[m],m=0...(2^d-1) are not permutation invariant.

The vertices should be ordered so that the
 macroelement (polyhedron) is mapped to the unit cube without
"twisting". This can be done by consistently labeling all faces. In
dimension d, if the vertices are numbered from 0 to (2^d-1) each
face corresponds to fixing exactly one bit in the binary
representation of the vertex label.

For example in 3D we have [000:111] as vertex labels, i.e. [0:7]
(this has nothing to do with coordinates) and the faces are:
(000,001,010,011)= fixing the first bit to be 0; and
(011,010,110,111)= fixing the second bit to be 1; and so on. The
ordering of the labels in the faces follows the same convention. In
the example we gave, the face for which the second bit is fixed is
numbered incorrectly, it should be (010,011,110,111).

POLAR COORDINATES have all angles in [0,pi] except the last one which
is in [0,2pi]

---

EOF README.md

---
