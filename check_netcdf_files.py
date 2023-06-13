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

try:
    ###################################################
    # Test only the last converted variable
    with open (sys.argv[-1], "r") as myfile:
        var=myfile.read().replace('\n', '').replace('copying variable', '').split(' ')
        test = test and len(var[6:])==len(nc.variables)
        print('testing', var[-1])
        test = test and np.isnan(nc.variables[var[-1]][:].sum()) or nc.variables[var[-1]][:].sum()==znc.variables[var[-1]][:].sum()
except:
    ###################################################
    # Test all variables:
    # print 'testing all vars'
    # for var in nc.variables.keys():
    #     test = test and np.isnan(nc.variables[var][:].sum()) or nc.variables[var][:].sum()==znc.variables[var][:].sum()
    test = False

###################################################
nc.close(); znc.close()
###################################################


if test: #files are identical
    os.remove(sys.argv[1])
    #os.rename(sys.argv[1],sys.argv[1]+'.trash')
    os.rename(sys.argv[2],sys.argv[1])

