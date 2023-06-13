#!/usr/bin/env python


##############################
from __future__ import print_function
import sys,os
from netCDF4 import Dataset
import numpy as np
##############################

#path, file = os.path.split(sys.argv[1])

nc = Dataset(sys.argv[1], 'r')
znc = Dataset(sys.argv[2], 'r')

#check all variables are in the destination files
test = len(znc.variables)==len(nc.variables)

number = [int, float, complex, np.int16, np.int32, np.int64,np.uint16, np.uint32, np.uint64, np.float16, np.float32, np.float64, np.complex64]

try:
    ###################################################
    # Test only the last converted variable
    with open (sys.argv[3], "r") as myfile:
        vars=myfile.read().replace('\n', '').replace('copying variable', '').replace('copying global attributes ..', '').replace('copying dimensions ..', '').split(' ')[1:]
        #all variables have been converted:
        test = test and len(vars)==len(nc.variables); 
        test1 = False
        #check sum for last converted one (if not all masked)
        for var in vars:
            if nc.variables[var].size>0 and nc.variables[var].dtype in number:
                try:
                    if not nc.variables[var][:].mask.all(): raise Exception('testing', var)
                except:
                    test1 = np.isnan(nc.variables[var][:].sum()) or nc.variables[var][:].sum()==znc.variables[var][:].sum()
                    break
        test = test and test1
except:
    ###################################################
    # Test all variables:
    print('testing all vars')
    for var in list(nc.variables.keys()):
        print('testing', var, 'size', nc.variables[var].size)
        if nc.variables[var].size>0 and nc.variables[var].dtype in number:
            try:
                if not nc.variables[var][:].mask.all(): raise Exception('testing', var)
            except:
                test = test and np.isnan(nc.variables[var][:].sum()) or nc.variables[var][:].sum()==znc.variables[var][:].sum()
        if not test: break
    #test = False

###################################################
nc.close(); znc.close()
###################################################


if test: #files are identical
    print('Files are identical __ replace .nc with .znc')
    os.rename(sys.argv[1],sys.argv[1]+'.trash')
    os.rename(sys.argv[2],sys.argv[1])
    os.remove(sys.argv[1]+'.trash')
else:
    print('Keep original .nc')
    os.remove(sys.argv[2])
