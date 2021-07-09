# ------------------ link swig files ---------------- #
from sys import path as psys
from getpass import getuser as gu
# fenics_hazmath_path = 'Dropbox/SIMULA/fenics_hazmath/'
fenics_hazmath_path = 'add/your/path/here/'
swig_hazmath_path = 'hazmath-1.0.0/swig_files'

the_dir = ''.join(['/home/', gu(),'/',fenics_hazmath_path,swig_hazmath_path])
print(the_dir)
# path to search for local modules appended
psys.append(the_dir)
# --------------------------------------------------- #
from haznics import *

# question: are strings in hazmath null char ended? 
# e.g. inifile is defined as char[256] 

def set_params(d, in_param): 
    for key in d: 
        if isinstance(d[key], str): 
            exec("in_param.%s = \"%s\"" % (key, d[key]))
        else: 
            exec("in_param.%s = %s" % (key, d[key]))


in_param = input_param() 

print (20*"-")
print ("default parameters created when just constructing the struct from Python, printing just a few") 
print (in_param.inifile)
print (len(in_param.inifile))
print (in_param.FE_type)
print (20*"-")

print ("setting the inifile parameter in Python ")
in_param.inifile = "input.dat"
print (in_param.inifile)
print (20*"-")


print (20*"-") 
print ("setting via a dict, kind of ugly via exec, but anyways ") 
d = {"print_level" : 1, "gridfile" : "mesh", "nquad" : 4, "trouble" : 3, "AMG_maxit" : 12  }

set_params(d, in_param) 

print (in_param.nquad) 
print (in_param.gridfile) 
print (in_param.trouble) 

amgparam = AMG_param()
param_amg_init(amgparam)
print (20*"-")
print("AMG param init: ")
param_amg_print(amgparam)
print (20*"-")

param_amg_set(amgparam, in_param)
print (20*"-") 
print ("AMG_param is able to read from in_param even though I have added the variable \"trouble\"  ") 
param_amg_print(amgparam)
print (20*"-")

print (" testing input of wrong type, AMG_maxit is now float " ) 
d = {"print_level" : 1, "gridfile" : "mesh", "nquad" : 4, "trouble" : 3, "AMG_maxit" : 12.0  }
set_params(d, in_param) 

