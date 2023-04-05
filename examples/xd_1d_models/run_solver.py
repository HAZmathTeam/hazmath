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

default_dim=2;

getopt = OptionParser(usage=usage)

getopt.add_option("-f", "--solver_input", action="store", type="string", dest="solver_input_file", help="solver input", metavar="SOLVER_INPUT",default="./input/solver.input")

getopt.add_option("-i", "--input_matrices_dir", action="store", type="string", dest="matrices_dir", help="Directory with matrices(npy)", metavar="MATRICES_DIR",default="./input/1d_matrices_"+str(default_dim)+"d/")

(op, args) = getopt.parse_args()

import ctypes
import os
op.solver_input_file=os.path.abspath(op.solver_input_file)
op.matrices_dir=os.path.abspath(op.matrices_dir)+'/'

libxd_1d_=ctypes.CDLL(os.path.abspath('./libxd_1d_.so'))

sfile=bytes(op.solver_input_file,'utf-8')
mdir=bytes(op.matrices_dir,'utf-8')
print("\nsolver input file: ",op.solver_input_file,"\nmatrices_dir=",op.matrices_dir,"\n")
libxd_1d_.solver_xd_1d(sfile,mdir)


