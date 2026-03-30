# AMR Grid Generation (HAZmath)

## Build and run

```sh
cd HAZMATH_DIR
make distclean; make config shared=yes suitesparse=yes; make install
cd examples/amr_grids
make clean; make
./amr_grids.ex input/3d_fichera.input
```

## Output files

Output filenames are generated automatically from the input filename:

    input/3d_fichera.input  +  num_refinements{3}  +  refinement_type{20}

produces:

    output/3d_fichera_rl3_rt20.msh   (Gmsh MSH 2.0; custom element types for dim >= 4
                                      because Gmsh MSH 2.0 only defines elements up to dim 3)
    output/3d_fichera_rl3_rt20.vtu   (VTK/ParaView, dim < 4 only)

The `.msh` file includes boundary faces (edges in 2D, triangles in 3D,
tetrahedra in 4D, etc.) with their boundary codes, followed by volume elements.

## Input file format

Rules:
1. Lines after `%` are comments; blank lines are ignored.
2. Variables: `variable{value}` or `variable{v1 v2 ...}` (space-separated arrays).

### Parameters

| Parameter | Description |
|---|---|
| `title{...}` | Title string |
| `dimension{d}` | Spatial dimension |
| `print_level{p}` | Output verbosity (0 = quiet) |
| `refinement_type{rt}` | See table below |
| `num_refinements{n}` | Number of refinement levels |
| `amr_marking_type{m}` | See marking types below |
| `num_refine_points{n}` | Number of points for marking types 33/34 |
| `data_refine_points{...}` | Point coordinates for marking types 33/34 (format: `csys type x y [z]`) |
| `err_stop_refinement{e}` | Stopping tolerance for AMR (not used in examples) |

### Refinement types (`refinement_type`)

| Value | Method | Description |
|---|---|---|
| 0 | DGS bisection | Adaptive bisection refinement (any dimension). Uses vertex coloring from Diening-Gehring-Storn (2025). |
| 20 | Uniform Bey | Freudenthal uniform refinement: each simplex is subdivided into 2^d children (any dimension). |
| 21 | Marked Bey (test) | Selective Bey refinement on odd-indexed simplices with conforming closure. Checkerboard test pattern. |
| 22 | Marked Bey (corner) | Selective Bey refinement on simplices near the origin, with conforming closure. Threshold shrinks per level for graded meshes (e.g., re-entrant corners). |

When `amr_marking_type{0}` and `refinement_type` is 0, DGS bisection
refines all simplices uniformly. When `refinement_type` is 20, Bey's
uniform refinement subdivides every simplex into 2^d children.

### Refinement algorithms

**Bey (Freudenthal) uniform refinement** (`refinement_type{20}`).
Midpoints are inserted on every edge. Each d-simplex is subdivided
into 2^d children using the Freudenthal partition described in [Bey 2000].
The children are enumerated by binary vectors s in {0,...,2^d-1}.
For a parent simplex with vertices v_0,...,v_d, child s has vertices
w_0,...,w_d where each w_k is either an original vertex v_i or the
midpoint m(v_i,v_j) of an edge, determined by a recurrence on
L_k, U_k driven by the bits of s. When applied to all simplices
simultaneously, the face patterns match across neighbors and the
resulting mesh is conforming. Works in any spatial dimension.

**DGS bisection** (`refinement_type{0}`).
Newest-vertex bisection using the generalized vertex coloring of
[Diening, Gehring, Storn 2025]. On the initial mesh, vertices are
colored with N+1 greedy colors and each simplex is reordered by
decreasing vertex generation number. Bisection splits a simplex
into two children across its tagged edge (determined by Algorithm 4
of [DGS 2025]). Conforming closure propagates bisections to
neighbors sharing the tagged edge. Works for any conforming initial
triangulation in any dimension, including adaptive refinement with
marking.

### Combining Bey and DGS: selective Bey with conforming closure

The selective Bey refinement (`refinement_type{21}` or `{22}`) combines
both algorithms. Only marked simplices are Bey-refined; the rest of
the mesh is closed to restore conformity using a two-phase procedure:

1. **Face-Bey closure.** A neighbor sharing a full (d-1)-face with a
   Bey-refined simplex has all C(d,2) edges of that face midpointed.
   Simple bisection cannot match the Bey face pattern (the 2^(d-1)
   sub-faces include a central simplex whose vertices are all
   midpoints, unreachable by bisection cuts which always connect a
   midpoint to an existing vertex). Instead, the (d-1)-dimensional
   Bey subdivision is applied to the non-conforming face and each
   sub-face is coned to the apex vertex, producing 2^(d-1) children
   whose faces match the Bey pattern exactly.

2. **Bisection closure.** After face-Bey closure, remaining
   non-conformity consists of single-edge hanging nodes (an edge with
   a midpoint not yet absorbed by the simplex). These are resolved by
   simple bisection: replace the simplex with two children, one for
   each half of the split edge. Bisection propagates along the strip
   of simplices sharing the affected edge, exactly as in standard
   DGS closure.

The loop iterates until no non-conforming simplices remain. The
algorithm works in any spatial dimension and has been verified on
2D through 5D meshes.

### Marking types (`amr_marking_type`)

| Value | Method | Description |
|---|---|---|
| 0 | Controlled by `refinement_type` | Uses `refinement_type` to select DGS (0), uniform Bey (20), or selective Bey (21/22) |
| 1 | Solve/estimate/mark/refine | Generic adaptive loop with user-defined solve, estimate, and mark functions |
| 33 | DGS near points | Mark simplices containing specified points, refine with DGS bisection |
| 34 | Bey+DGS near points | Mark simplices near specified points (barycenter distance with shrinking threshold), refine with selective Bey + face-Bey/bisection closure |
| 35 | Bey then DGS near points | Bey refinement near points, then DGS adaptive completion |
| 44 | Features from file | Refine around features read from an external data file |

Types 33, 34, and 35 use `num_refine_points` and `data_refine_points` from
the input file. Types 34 and 35 start with threshold 1.0 and shrink per
level, producing graded meshes concentrated near the specified points.

**Example** — Fichera corner refined near the re-entrant point (0,0,0):
```
amr_marking_type{34}
num_refinements{6}
num_refine_points{1}
data_refine_points{0 0   0.0 0.0   0.0}
```

### Geometry specification

**Coordinate systems** (`num_coordsystems`, `data_coordsystems`):
Each row = `label  type  origin[0:dim-1]`. Type 0 = Cartesian, 1 = polar.

**Vertices** (`num_vertices`, `data_vertices`):
Each row = `label  coord_system  x[0:dim-1]`.

**Edges** (`num_edges`, `data_edges`):
Each row = `v1  v2  num_divisions`. Parallel edges (mapped to the unit cube)
share the maximum division count. Missing edges default to 1 division.

**Macroelements** (`num_macroelements`, `data_macroelements`):
Each row = `v[0] ... v[2^d-1]  material_code`. Vertices must be ordered so
the macroelement maps to the unit cube without singularity (see below).

**Macrofaces** (`num_macrofaces`, `data_macrofaces`):
Each row = `v[0] ... v[2^(d-1)-1]  boundary_code`. A macroface can be
on the geometric boundary or shared between two macroelements (interior
face with a prescribed code, e.g., material interface). Unlisted
boundary faces get a default code `(face_number % 8) + 1`; unlisted
interior faces get code 0.

### Vertex ordering for macroelements

Vertices `v[0]...v[2^d-1]` are mapped to the unit cube vertices
`0...2^d-1` (binary coordinates). The ordering is consistent if every
k-dimensional face maps to a k-dimensional cube face (vertices sharing
d-k fixed bits). The resulting d-linear map must have nonzero Jacobian.

## Input files included

| File | Title |
|------|-------|
| `input/2d_2L.input` | 2 Lshaped domains with connected boundaries |
| `input/2d_ann.input` | Annulus in 2D |
| `input/2d_circle.input` | Circle in 2D |
| `input/2d_grid.input` | 2D Sector (polar coords) |
| `input/2d_L.input` | L-shaped domain |
| `input/2d_L_adaptive.input` | L-shaped domain (adaptive near corner) |
| `input/2d_SQ+L.input` | Grid with two connected components and connected bndry |
| `input/2d_square.input` | Unit square in 2D |
| `input/3d_2cubes_edge.input` | Two cubes sharing an edge |
| `input/3d_2cubes_vertex.input` | Two cubes sharing a vertex |
| `input/3d_cube.input` | Grid on the cube (0,1)^3 |
| `input/3d_fichera.input` | Fichera corner |
| `input/3d_fichera_adaptive.input` | Fichera corner (adaptive near corner) |
| `input/3d_fichera_bey_dgs.input` | Fichera corner (Bey + DGS) |
| `input/4d_cube.input` | Grid on the 4d cube (-1,1)^4 |
| `input/5d_cube.input` | Grid on the 5d cube (-1,1)^5 |

## Source files

| File | Contents |
|---|---|
| `src/amr/uniform_refinement.c` | `get_edge_nd`, `uniformrefine` (Bey/Freudenthal, any dim), `uniformrefine_marked` (selective Bey + closure) |
| `src/amr/amr_core.c` | `refine` (DGS bisection), `haz_refine_simplex`, `haz_bisect_new/reuse`, `make_uniform_mesh` |
| `src/amr/scomplex.c` | Simplicial complex: init, free, geometry, volumes, FEM data, boundary, conformity check |

## Face codes and boundary data

After `sc_build_fem_data(sc)`, face codes are stored in `sc->fem->f_flag[f]`
for every face in the mesh. Faces with nonzero codes are listed in
`sc->fem->coded_faces[]` with a boundary/interior indicator in
`sc->fem->coded_f_btype[]` (0 = boundary, 1 = interior).

Face codes are inherited from the macroface definitions in the input
file. A mesh face has code C if all its vertices belong to macroface C
(tracked via `sc->bndry_v`). Interior faces between macroelements can
have nonzero codes (e.g., material interfaces, re-entrant cavity walls).
All coded faces (boundary + interior) survive Gmsh `.msh` round-trips.

Boundary faces without a prescribed code get a default of
`(face_number % 8) + 1`, ensuring every boundary face has a nonzero code.

Vertex codes (`sc->bndry[v]`) are the minimum of the face codes meeting
at vertex `v`, computed in `scfinalize()`.

### Boundary condition convention (`include/fem.h`)

| Range | Macro | Meaning |
|-------|-------|---------|
| 0 | — | Interior (no BC) |
| 1-16 | `MARKER_DIRICHLET` | Dirichlet |
| 17-32 | `MARKER_NEUMANN` | Neumann |
| 33-64 | `MARKER_ROBIN` | Robin |
| 65+ | `MARKER_BOUNDARY_NO` | No BC applied |

See `HAZMATH_MESH_RW.md` for C code examples of looping over coded
faces and applying Dirichlet, Neumann, and Robin conditions using
`simplex_local_data` and `fe_local_data`.

The test program `test_bndry_codes.c` demonstrates writing boundary
faces as VTU with face codes for ParaView visualization.

## References

[Bey 2000] J. Bey, *Simplicial grid refinement: on Freudenthal's
algorithm and the optimal number of congruence classes*,
Numer. Math., 85(1):1-29, 2000.
DOI: [10.1007/s002110050477](https://doi.org/10.1007/s002110050477)

[DGS 2025] L. Diening, L. Gehring, J. Storn, *Adaptive Mesh Refinement
for Arbitrary Initial Triangulations*, Found. Comput. Math., 2025.
DOI: [10.1007/s10208-024-09642-1](https://doi.org/10.1007/s10208-024-09642-1)

[Freudenthal 1942] H. Freudenthal, *Simplizialzerlegungen von
beschrankter Flachheit*, Ann. of Math., 43(3):580-582, 1942.
