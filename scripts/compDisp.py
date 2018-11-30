#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 13:52:36 2018

@author: guru

Normalize kappa0 and kappaL data: for the curves reported in M. Schrade, see if
kappa_L over kappa_0 is the same for each curve

Debye curve was hard to data thief... kappa0 is probably incorrect
"""

import numpy as np

#Halfnium sample

fileFit = "/Users/guru/Documents/PDscattering/datafiles/FitDisp_Hf.csv"
fileDebye = "/Users/guru/Documents/PDscattering/datafiles/Debye_Hf.csv"
fileDebye2 = "/Users/guru/Documents/PDscattering/datafiles/Debye_max_Hf.csv"

dFit = np.genfromtxt(fileFit, dtype = 'float', delimiter = ',')
cFit = dFit[:,0]
kFit = dFit[:,1]
k0Fit = kFit[0]

dDebye = np.genfromtxt(fileDebye, dtype = 'float', delimiter = ',')
cDebye = dDebye[:,0]
kDebye = dDebye[:,1]
k0Debye = kDebye[0]

dDebye2 = np.genfromtxt(fileDebye2, dtype = 'float', delimiter = ',')
cDebye2 = dDebye2[:,0]
kDebye2 = dDebye2[:,1]
k0Debye2 = kDebye2[0]


print(kFit/k0Fit)
print(kDebye/k0Debye)
print(kDebye2/k0Debye2)


#Zirconium sample


fileFit = "/Users/guru/Documents/PDscattering/datafiles/FitDisp_Zr.csv"
fileDebye = "/Users/guru/Documents/PDscattering/datafiles/Debye_Zr.csv"
fileDebye2 = "/Users/guru/Documents/PDscattering/datafiles/Debye_max_Zr.csv"

dFit = np.genfromtxt(fileFit, dtype = 'float', delimiter = ',')
cFit = dFit[:,0]
kFit = dFit[:,1]
k0Fit = kFit[0]

dDebye = np.genfromtxt(fileDebye, dtype = 'float', delimiter = ',')
cDebye = dDebye[:,0]
kDebye = dDebye[:,1]
k0Debye = kDebye[0]

dDebye2 = np.genfromtxt(fileDebye2, dtype = 'float', delimiter = ',')
cDebye2 = dDebye2[:,0]
kDebye2 = dDebye2[:,1]
k0Debye2 = kDebye2[0]


print(kFit/k0Fit)
print(kDebye/k0Debye)
print(kDebye2/k0Debye2)
