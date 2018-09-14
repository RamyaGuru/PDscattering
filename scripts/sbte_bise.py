#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 17:08:56 2018

@author: guru

Bi2Se3-Sb2Te3 test file

Plot the gammas and see how much they differ for both examples of the monatomic
lattice approximation
"""
from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/guru/Documents/PDscattering/datafiles/"

#Compute the thermal conductivity for both cases

prop = PropsM(directory + "mass_cancel.csv")

crange = np.linspace(0.00001,0.9999,5)

gammaS = np.zeros(5)
gammaUC = np.zeros(5)
gammaUC2 = np.zeros(5)
i=0
for c in crange:
    prop.bfz(c)
    gammaS[i] = prop.gamma
    prop.puc(c)
    gammaUC[i] = prop.gamma
    prop.puc_proper2(c)
    gammaUC2[i] = prop.gamma 
    i = i+1

plt.plot(crange, gammaS, label = 'bfz')
plt.plot(crange, gammaUC, label = 'puc')
plt.plot(crange, gammaUC2, label = 'puc2')
plt.legend()

