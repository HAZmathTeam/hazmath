
# INPUT FILES for the simple grid generator in HAZMATH

---

HAZMATH: A Simple Finite Element, Graph, and Solver Library

Copyright (c) 2009- HAZMath: Xiaozhe Hu, James H Adler, Ludmil T Zikatanov  

---

Here we construct 2d and 3d meshes around 1d graphs. A python script calls the hazmath refinement routine

The usage is: python3 run_meshes_xd_1d.py [-h] -d dimension -i input_dir -o output_dir 

The defaults are: dim=3; input dir is ./input/1d_nets_Xd for dimension X (2 or 3) and output dir is output.

Try "python3 run_meshes_xd_1d.py -h" for more on usage;

The filenames in input/* are fixed and cannot be changed unless you modify the C-source. But there is an option to take any file from anywhere else in place of these. Indeed,

Try "python3 run_meshes_xd_1d.py -h" for more on usage;


---

EOF README.md

---
