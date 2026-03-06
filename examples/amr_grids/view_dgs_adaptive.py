"""
view_dgs_adaptive.py

Visualize 2D and 3D DGS adaptive meshes using matplotlib.
Reads VTU files and overlays refinement points.

Usage:
  python3 view_dgs_adaptive.py
  (run test_dgs_adaptive.ex first to generate output/*.vtu and output/*.csv)
"""
import os, csv, struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import xml.etree.ElementTree as ET
import base64

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

def read_vtu(filename):
    """Parse a VTK XML UnstructuredGrid (.vtu) file. Returns points, cells."""
    tree = ET.parse(filename)
    root = tree.getroot()
    piece = root.find('.//Piece')
    npts = int(piece.attrib['NumberOfPoints'])
    ncells = int(piece.attrib['NumberOfCells'])

    # Find Points DataArray
    pts_da = piece.find('.//Points/DataArray')
    fmt = pts_da.attrib.get('format', 'ascii')
    if fmt == 'ascii':
        pts_text = pts_da.text.strip().split()
        pts_raw = np.array([float(x) for x in pts_text])
    else:
        # binary or appended — fallback
        pts_raw = np.fromstring(base64.b64decode(pts_da.text.strip()), dtype=np.float64)
    ncomp = int(pts_da.attrib.get('NumberOfComponents', 3))
    points = pts_raw.reshape(-1, ncomp)

    # Find Cells DataArrays
    cells_section = piece.find('.//Cells')
    conn_da = None
    off_da = None
    for da in cells_section.findall('DataArray'):
        if da.attrib['Name'] == 'connectivity':
            conn_da = da
        elif da.attrib['Name'] == 'offsets':
            off_da = da

    conn_text = conn_da.text.strip().split()
    connectivity = np.array([int(x) for x in conn_text])
    off_text = off_da.text.strip().split()
    offsets = np.array([int(x) for x in off_text])

    # Build cell list
    cells = []
    prev = 0
    for o in offsets:
        cells.append(connectivity[prev:o].tolist())
        prev = o

    return points, cells

def read_points_csv(csv_file):
    pts = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            pts.append([float(x) for x in row])
    return pts

def visualize_2d(vtu_file, csv_file, png_file):
    print("Visualizing 2D: %s" % vtu_file)
    points, cells = read_vtu(vtu_file)
    x = points[:, 0]
    y = points[:, 1]

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    # Draw triangles
    triangles = np.array([c for c in cells if len(c) == 3])
    triang = mtri.Triangulation(x, y, triangles)
    ax.triplot(triang, 'k-', linewidth=0.3)
    ax.tripcolor(triang, np.zeros(len(x)), cmap='Blues', alpha=0.15,
                 edgecolors='black', linewidth=0.3)

    # Draw refinement points
    if os.path.exists(csv_file):
        pts = read_points_csv(csv_file)
        colors = ['red', 'green', 'blue']
        for i, p in enumerate(pts):
            ax.plot(p[0], p[1], 'o', color=colors[i % 3], markersize=8,
                    markeredgecolor='black', markeredgewidth=0.5, zorder=10)
            ax.annotate('P%d' % i, (p[0]+0.012, p[1]+0.012), fontsize=9,
                        color=colors[i % 3], fontweight='bold')

    ax.set_aspect('equal')
    ax.set_title('DGS Adaptive Refinement (2D)\n%d triangles, %d vertices' %
                 (len(triangles), len(x)), fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.tight_layout()
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % png_file)

def visualize_3d(vtu_file, csv_file, png_file):
    print("Visualizing 3D: %s" % vtu_file)
    points, cells = read_vtu(vtu_file)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Collect boundary faces: faces that appear only once
    face_count = {}
    for cell in cells:
        if len(cell) != 4:
            continue
        # 4 faces of a tetrahedron
        faces = [
            tuple(sorted([cell[1], cell[2], cell[3]])),
            tuple(sorted([cell[0], cell[2], cell[3]])),
            tuple(sorted([cell[0], cell[1], cell[3]])),
            tuple(sorted([cell[0], cell[1], cell[2]])),
        ]
        for f in faces:
            face_count[f] = face_count.get(f, 0) + 1

    bnd_faces = [f for f, c in face_count.items() if c == 1]
    print("  Boundary faces: %d" % len(bnd_faces))

    # Draw boundary faces
    polys = []
    for f in bnd_faces:
        verts = [points[i] for i in f]
        polys.append(verts)

    pc = Poly3DCollection(polys, alpha=0.25, facecolor='lightskyblue',
                          edgecolor='gray', linewidth=0.3)
    ax.add_collection3d(pc)

    # Draw refinement points
    if os.path.exists(csv_file):
        pts = read_points_csv(csv_file)
        colors = ['red', 'green', 'blue']
        for i, p in enumerate(pts):
            ax.scatter(p[0], p[1], p[2], color=colors[i % 3], s=80,
                       edgecolors='black', linewidths=0.5, zorder=10,
                       depthshade=False)

    ax.set_xlim(x.min()-0.05, x.max()+0.05)
    ax.set_ylim(y.min()-0.05, y.max()+0.05)
    ax.set_zlim(z.min()-0.05, z.max()+0.05)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('DGS Adaptive Refinement (3D)\n%d tetrahedra, %d vertices' %
                 (len(cells), len(x)), fontsize=14)
    ax.view_init(elev=20, azim=35)
    plt.tight_layout()
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % png_file)

# Run
vtu2 = os.path.join(outdir, "dgs_adaptive_2d.vtu")
csv2 = os.path.join(outdir, "dgs_points_2d.csv")
png2 = os.path.join(outdir, "dgs_adaptive_2d.png")
if os.path.exists(vtu2):
    visualize_2d(vtu2, csv2, png2)

vtu3 = os.path.join(outdir, "dgs_adaptive_3d.vtu")
csv3 = os.path.join(outdir, "dgs_points_3d.csv")
png3 = os.path.join(outdir, "dgs_adaptive_3d.png")
if os.path.exists(vtu3):
    visualize_3d(vtu3, csv3, png3)

print("\nDone.")
