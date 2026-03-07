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

    output/3d_fichera_rl3_rt20.haz   (HAZmath grid format)
    output/3d_fichera_rl3_rt20.msh   (Gmsh MSH 2.0, dim < 4 only)
    output/3d_fichera_rl3_rt20.vtu   (VTK/ParaView, dim < 4 only)

The `.msh` file includes boundary faces (edges in 2D, triangles in 3D)
with their boundary codes, followed by volume elements.

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
| `refinement_type{rt}` | 20 = DGS + criss-cross; 21 = DGS + consistent diagonals |
| `num_refinements{n}` | Number of refinement levels |
| `amr_marking_type{m}` | 0 = uniform; nonzero = user-defined marking |
| `err_stop_refinement{e}` | Stopping tolerance for AMR (not used in examples) |

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
Each row = `v[0] ... v[2^(d-1)-1]  boundary_code`. Unlisted boundary faces
get code 1 (Dirichlet); unlisted interior faces get code 0.

### Vertex ordering for macroelements

Vertices `v[0]...v[2^d-1]` are mapped to the unit cube vertices
`0...2^d-1` (binary coordinates). The ordering is consistent if every
k-dimensional face maps to a k-dimensional cube face (vertices sharing
d-k fixed bits). The resulting d-linear map must have nonzero Jacobian.

## Input files included

```
input/2d_2L.input           input/3d_2cubes_edge.input
input/2d_ann.input          input/3d_2cubes_vertex.input
input/2d_circle.input       input/3d_cube.input
input/2d_grid.input         input/3d_fichera.input
input/2d_L.input            input/4d_cube.input
input/2d_SQ+L.input         input/5d_cube.input
input/2d_square.input
```

---
Previous README: [README.md.bak](README.md.bak)
