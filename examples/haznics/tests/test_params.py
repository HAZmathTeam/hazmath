from haznics import *

print("Input parameters with constructor that sets default by calling init function")
in_param = input_param()

print("in_param type", in_param)
print("in_param dir", dir(in_param))
print("in_param vars", vars(in_param))

print(20*"-")
print("Default parameters created when just constructing the struct from Python, printing just a few")
print(in_param.inifile)
print(len(in_param.inifile))
print(in_param.FE_type)
print(in_param.AMG_smoother)
print(in_param.linear_itsolver_type)
print(20*"-")

print("Input parameters with constructor that sets default from input file")
in_param2 = input_param("input.dat")

print(20*"-")
print("Default parameters created when just constructing the struct from Python, printing just a few:")
print(in_param2.inifile)
print(len(in_param2.inifile))
print(in_param2.FE_type)
print(in_param2.AMG_smoother)
print(in_param2.linear_itsolver_type)
print(20*"-")

print("Setting the init file parameter in Python ")
in_param.inifile = "input.dat"
print(in_param.inifile)
print(20*"-")

print(20*"-")
print("Setting input parameters via a dict")
d = {"print_level" : 1, "gridfile" : "mesh", "nquad" : 4, "trouble" : 3, "AMG_maxit" : 12  }
param_input_set_dict(d, in_param)

print("Checking some input parameters: ")
print(in_param.nquad)
print(in_param.gridfile)
print(in_param.trouble)

print("AMG params with constructor that sets default by calling init function")
amgparam = AMG_param()
print(20*"-")
print("AMG param after init: ")
param_amg_print(amgparam)
print(20*"-")

print("Setting AMG parameters from input parameters")
param_amg_set(amgparam, in_param)
print(20*"-")
print("AMG_param is able to read from in_param even though I have added the variable \"trouble\"  ")
param_amg_print(amgparam)
print(20*"-")

print(" Testing input of wrong type, AMG_maxit is now float " )
d = {"print_level" : 1, "gridfile" : "mesh", "nquad" : 4, "trouble" : 3, "AMG_maxit" : 12.0  }
param_input_set_dict(d, in_param)

