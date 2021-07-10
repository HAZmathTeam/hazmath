from sys import path as psys
from os import path as pos
# (ASSUMING that we  are in HAZMATH_DIR/examples/haznics):
# path to search for local modules 
the_dir = ''.join([pos.realpath('../../'),'/swig_files/'])
print("Adding to syspath:",the_dir);
psys.append(the_dir)
