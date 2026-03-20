# HAZMATH Mesh Read/Write

## Reading and Writing Meshes

HAZMATH stores meshes as simplicial complexes (`scomplex`). The
primary mesh format is Gmsh `.msh` v2 ASCII. VTK/VTU is supported
for visualization output. MATLAB `.m` output is available for quick
plotting.

### Reading a Mesh

```c
#include "hazmath.h"

// Read a Gmsh .msh file
scomplex *sc = sc_read_gmsh("mesh.msh");

// Build mesh connectivity and FEM data
find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
sc_vols(sc);
sc_build_fem_data(sc);
```

`sc_read_gmsh()` returns a `scomplex*` with:
- `sc->nodes` — element-to-vertex connectivity
- `sc->x` — vertex coordinates (row-major: `x[i*dim+d]`)
- `sc->flags` — element material/region tags
- `sc->bndry` — per-vertex boundary codes (min over adjacent faces)
- `sc->bndry_f2v` — boundary face-to-vertex map with face codes
- `sc->bndry_v` — vertex-to-boundary-face map with face codes

After reading, call `find_nbr`, `sc_vols`, and `sc_build_fem_data` to
populate neighbor lists, element volumes, and FEM data (edges, faces,
quadrature, etc.).

### Writing a Mesh

```c
// Write to Gmsh .msh format
sc_write_gmsh("output.msh", sc, 1);

// Write to VTK/VTU format (for ParaView, VisIt, etc.)
vtu_data vdata;
vtu_data_init(sc, &vdata);
sc_write_vtk("output.vtu", &vdata);
vtu_data_free(&vdata);

// Write to MATLAB .m format (for quick plotting)
sc_matlab_write(sc, "output.m");
```

### Function Reference

| Function | File | Description |
|----------|------|-------------|
| `sc_read_gmsh(filename)` | `src/utilities/io.c` | Read Gmsh .msh v2 ASCII, return `scomplex*` |
| `sc_write_gmsh(filename, sc, shift)` | `src/utilities/io.c` | Write `scomplex` to Gmsh .msh v2 ASCII |
| `sc_write_vtk(filename, vdata)` | `src/amr/amr_utils.c` | Write mesh + data to VTK/VTU format |
| `sc_matlab_write(sc, filename)` | `src/utilities/io.c` | Write mesh to MATLAB .m file |
| `vtu_data_init(sc, vdata)` | `src/amr/amr_utils.c` | Initialize VTU data with boundary/element codes |
| `vtu_data_free(vdata)` | `src/amr/amr_utils.c` | Free VTU data arrays |

### Dimension Support

`sc_read_gmsh` and `sc_write_gmsh` work for **any spatial dimension**,
including d > 3. For dimensions 1, 2, and 3 the standard Gmsh element
types are used (line, triangle, tetrahedron). For d > 3, a custom
element type `50 + d` is used to represent (d+1)-node simplices. This
allows round-tripping of higher-dimensional simplicial complexes
through the .msh format.

`sc_matlab_write` also works for any dimension.

`sc_write_vtk` is limited to dimensions 1, 2, and 3 (VTK/VTU format
does not support higher-dimensional elements).

### Generating Meshes Internally

```c
// Unit cube mesh in dimension dim, with refinement
scomplex *sc = make_uniform_mesh(dim, ref_levels, ref_type, set_bndry);

// 1D mesh on [a, b]
scomplex *sc;
create1Dgrid_Line(&sc, left, right, nelm);

// Mesh from domain graph description (any dimension)
scomplex **sc_all = mesh_cube_init(dim, ndiv, ref_type);
scomplex *sc = sc_all[0];
scfinalize(sc, NULL, 1);  // NULL = compact in-place, discard hierarchy
sc_vols(sc);
sc_build_fem_data(sc);
free(sc_all);

// To keep the hierarchy and extract leaf mesh separately:
// scomplex sc_leaf;
// scfinalize(sc, &sc_leaf, 1);  // sc keeps hierarchy, sc_leaf gets leaves
```

---

## Round-Trip Example: Two Cubes with Boundary Codes

The example `examples/amr_grids/test_msh_roundtrip.c` demonstrates the
full mesh I/O workflow with boundary code preservation through refinement.

**What it does:**
1. Generates a 3D mesh of two cubes sharing an edge from
   `input/3d_2cubes_edge.input` (12 distinct boundary face codes)
2. Writes to `output/2cubes_original.msh` and `.vtu`
3. Reads back from the `.msh` file using `sc_read_gmsh`
4. Verifies conformity of the imported mesh
5. Refines once using newest vertex bisection
6. Propagates boundary codes to new vertices via `find_cc_bndry_cc`
7. Writes refined mesh to `output/2cubes_refined.msh`, `.vtu`, and `.m`

**Running:**
```
cd examples/amr_grids
make
./test_msh_roundtrip.ex
```

**Boundary codes in the input file:** each boundary face of the two-cube
domain has a distinct code (1 through 12), allowing verification that
codes survive the write → read → refine cycle. After refinement, all
boundary vertices have nonzero codes derived from the face-level
boundary data via `sc->bndry_v`.

---

## Gmsh .msh v2 ASCII Format

This section describes the subset of the Gmsh format used by HAZMATH.
For the full specification, see:
https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

### File Structure

```
$MeshFormat
2.0 0 8
$EndMeshFormat
$Nodes
...
$EndNodes
$Elements
...
$EndElements
```

### Nodes Section

```
$Nodes
<number_of_nodes>
<node_id> <x> <y> [<z>]
...
$EndNodes
```

- `node_id` is a 1-based integer identifier
- Coordinates are floating-point; 2D meshes may omit `z` or set it to 0
- HAZMATH subtracts the minimum node ID to convert to 0-based indexing

**Example (4 nodes in 2D):**
```
$Nodes
4
1 0.0 0.0
2 1.0 0.0
3 1.0 1.0
4 0.0 1.0
$EndNodes
```

### Elements Section

```
$Elements
<total_number_of_elements>
<elm_id> <elm_type> <num_tags> <tag1> [<tag2> ...] <node1> <node2> ...
...
$EndElements
```

Each element has:
- `elm_id` — 1-based element identifier
- `elm_type` — integer specifying the element shape (see below)
- `num_tags` — number of integer tags following
- Tags: the first tag is typically the physical group (used as boundary
  code or material flag); the second tag (if present) is the elementary
  entity
- Node list: 1-based node IDs defining the element connectivity

### Element Types Used by HAZMATH

| Type | Shape | Nodes | Dimension |
|------|-------|-------|-----------|
| 15 | Point | 1 | 0 (boundary vertex in 1D) |
| 1 | Line | 2 | 1 (edge / boundary face in 2D) |
| 2 | Triangle | 3 | 2 (face / boundary face in 3D) |
| 4 | Tetrahedron | 4 | 3 (volume element) |
| 50+d | d-simplex | d+1 | d (HAZMATH custom, for d > 3) |

### Volume Elements vs Boundary Faces

A .msh file typically contains two kinds of elements:

1. **Volume elements** — the simplices of the mesh (lines in 1D,
   triangles in 2D, tetrahedra in 3D)
2. **Boundary face elements** — lower-dimensional elements on the
   boundary (points in 1D, lines in 2D, triangles in 3D)

HAZMATH determines the spatial dimension from the highest-dimensional
element type present. Volume elements become `sc->nodes` and
`sc->flags`. Boundary face elements become `sc->bndry_f2v` with face
codes from the first tag.

**Example (2 triangles with 1 boundary edge):**
```
$Elements
3
1 1 2 1 1 1 2
2 2 1 1 1 2 3
3 2 1 1 1 3 4
$EndElements
```

#### Field-by-field breakdown of `1 1 2 1 1 1 2`:

```
1  1  2  1  1  1  2
│  │  │  │  │  │  └─ node 2 (vertex 2, 1-based)
│  │  │  │  │  └──── node 1 (vertex 1, 1-based)
│  │  │  │  └─────── tag 2: elementary entity ID (geometric entity in Gmsh)
│  │  │  └────────── tag 1: physical group = 1 (boundary code used by HAZMATH)
│  │  └───────────── num_tags = 2 (two integer tags follow)
│  └──────────────── elm_type = 1 (2-node line segment)
└─────────────────── elm_id = 1 (1-based element identifier)
```

This defines a boundary edge between vertices 1 and 2 with **boundary
code 1** (the physical group tag). HAZMATH uses the first tag as the
boundary code for boundary face elements and as the material/region
flag for volume elements.

#### Field-by-field breakdown of `2 2 1 1 1 2 3`:

```
2  2  1  1  1  2  3
│  │  │  │  │  │  └─ node 3 (vertex 3)
│  │  │  │  │  └──── node 2 (vertex 2)
│  │  │  │  └─────── node 1 (vertex 1)
│  │  │  └────────── tag 1: physical group = 1 (material/region flag)
│  │  └───────────── num_tags = 1 (one tag)
│  └──────────────── elm_type = 2 (3-node triangle)
└─────────────────── elm_id = 2
```

This defines a triangle (volume element in 2D) with vertices 1, 2, 3
and material flag 1.

#### How HAZMATH distinguishes boundary faces from volume elements

The spatial dimension is determined from the highest-dimensional element
type present. Elements of that dimension become volume elements
(`sc->nodes`, `sc->flags`). Elements of dimension `dim-1` become
boundary faces (`sc->bndry_f2v`) with boundary codes from the first tag.
Lower-dimensional elements (e.g., points in a 3D mesh) are ignored.

---

## How Boundary Codes Work

### In the .msh File

Boundary faces carry a physical group tag (the first tag). For example,
to mark different parts of a square boundary:

Consider the unit square with vertices:
```
$Nodes
4
1 0.0 0.0       % vertex 1: bottom-left
2 1.0 0.0       % vertex 2: bottom-right
3 1.0 1.0       % vertex 3: top-right
4 0.0 1.0       % vertex 4: top-left
$EndNodes
```

The four boundary edges with distinct codes:
```
1 1 2 1 1 1 2    % edge: type=line, 2 tags, phys=1 elem=1, verts 1-2 → bottom, code 1
2 1 2 2 2 2 3    % edge: type=line, 2 tags, phys=2 elem=2, verts 2-3 → right,  code 2
3 1 2 3 3 3 4    % edge: type=line, 2 tags, phys=3 elem=3, verts 3-4 → top,    code 3
4 1 2 4 4 4 1    % edge: type=line, 2 tags, phys=4 elem=4, verts 4-1 → left,   code 4
```

Each line follows the pattern: `elm_id  elm_type  num_tags  phys_tag  elem_tag  node1  node2`.
The physical tag (first tag) becomes the HAZMATH boundary code.
HAZMATH uses boundary codes 1-16 for Dirichlet and 17-32 for Neumann by convention.

### In HAZMATH After Reading

`sc_read_gmsh()` stores boundary data in three structures:

1. `sc->bndry_f2v` — iCSRmat, rows = boundary faces, `JA` = vertex
   indices, `val` = face code (same for all vertices of a face)

2. `sc->bndry_v` — iCSRmat (transpose of `bndry_f2v`), rows = vertices,
   `JA` = boundary face indices, `val` = face codes. A corner vertex
   belonging to two boundary faces will have two entries.

3. `sc->bndry[v]` — per-vertex code, derived as the **minimum** code
   over all boundary faces touching vertex `v`. Interior vertices have
   code 0.

### After Refinement

When the mesh is refined, new boundary vertices inherit codes from their
parent vertices via `sc->bndry_v`. The function `find_cc_bndry_cc()`
computes the intersection of boundary face memberships of the two parent
vertices. This ensures boundary codes are preserved correctly through
any number of refinement levels.

---

## Generating .msh Files with Gmsh

To create a mesh with boundary codes in Gmsh:

1. Define physical groups for boundary curves/surfaces in the `.geo` file:
   ```
   Physical Line("bottom") = {1};
   Physical Line("right") = {2};
   Physical Surface("domain") = {1};
   ```

2. Generate the mesh and save as .msh v2:
   ```
   gmsh -2 -format msh2 -o mesh.msh geometry.geo
   ```

3. Read in HAZMATH:
   ```c
   scomplex *sc = sc_read_gmsh("mesh.msh");
   find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
   sc_vols(sc);
   sc_build_fem_data(sc);
   ```

---

## Reference

- Gmsh file format: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
- Gmsh software: https://gmsh.info
