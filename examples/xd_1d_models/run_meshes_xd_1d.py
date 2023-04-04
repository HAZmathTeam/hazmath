# -*- coding: utf-8 -*-
"""
Created at 16:25:34 on Sun 20230404

@author: ltz1
"""

# 
from optparse import OptionParser
import time
import re
import pathlib
import meshio

#ztime=str(time.strftime("%H:%M:%S"))
#zdate=str(time.strftime("%Y%m%d"))

##############################################################

usage = "usage: python3 %prog [-h] -d dimension -i input_dir -o output_dir"

default_dim=3;

getopt = OptionParser(usage=usage)

getopt.add_option("-d", "--dim", action="store", type="int",dest="dim", help="X in Xd-1d (default="+str(default_dim)+")", metavar="DIMENSION",default=default_dim)

getopt.add_option("-m", "--max_nodes", action="store", type="int",dest="max_nodes", help="max refinement levels"+str(default_dim)+")", metavar="MAX_NODES",default=10)

getopt.add_option("-r", "--ref_levels", action="store", type="int",dest="ref_levels", help="max refinement levels"+str(default_dim)+")", metavar="REF_LEVELS",default=100)

getopt.add_option("-i", "--input_dir", action="store", type="string", dest="idir", \
                  help="Input directory (see README.md)",metavar="INPUT_DIR",default="./input/1d_nets_"+str(default_dim)+"d/")

getopt.add_option("-o", "--output_dir", action="store", type="string", dest="odir", help="Output directory", metavar="OUTPUT_DIR",default="./output/")

(op, args) = getopt.parse_args()

if(op.dim<2):
    op.dim=2
elif(op.dim>3):
    op.dim=3

import ctypes
import os
op.idir=os.path.abspath(op.idir)+'/'
op.odir=os.path.abspath(op.odir)+'/'

sdim=str(op.dim)
libxd_1d_=ctypes.CDLL(os.path.abspath('./libxd_1d_.so'))

print("\n%%%%%%DIM=",op.dim,\
        "\nmax_nodes=",op.max_nodes, \
      "\nref_levels=",op.ref_levels,\
      "\ninput_dir=",op.idir,\
	"\noutput_dir=",op.odir,"\n");

dim=ctypes.c_int(op.dim);
max_nodes=ctypes.c_int(op.max_nodes);
ref_levels=ctypes.c_int(op.ref_levels);

#idir=ctypes.c_wchar_p(op.idir);
#odir=ctypes.c_wchar_p(op.odir);
idir=bytes(op.idir,'utf-8')
odir=bytes(op.odir,'utf-8')


#libxd_1d_.xd_1d_lib(const INT dimbig, const INT max_nodes_in, const INT ref_levels_in, const char *idir, const char *odir)
libxd_1d_.xd_1d_lib(dim,max_nodes,ref_levels,idir,odir)

mesh=meshio.read(op.odir+'1d_grid.vtu')
mesh.write(op.odir+'1d_grid.xdmf')

mesh=meshio.read(op.odir+sdim+'d_grid.vtu')
mesh.write(op.odir+sdim+'d_grid.xdmf')




#https://github.com/MiroK/metric-amg-examples/blob/master/src/emi3d1d_sylvie_ludmil.py


####################################################################
# LATER
####################################################################
# mesh3d = Mesh()
#     with HDF5File(mesh3d.mpi_comm(), './data/sylvie/sylvie_3d.h5', 'r') as h5:
#         h5.read(mesh3d, '/mesh', False)

#     mesh1d = Mesh()
#     with HDF5File(mesh1d.mpi_comm(), './data/sylvie/sylvie_1d.h5', 'r') as h5:
#         h5.read(mesh1d, '/mesh', False)

#     mesh1d_radii = MeshFunction('double', mesh1d, 1, 0)
#     with HDF5File(mesh1d.mpi_comm(), './data/sylvie/sylvie_1d.h5', 'r') as h5:    
#         h5.read(mesh1d_radii, '/radii_cell_foo')

        
