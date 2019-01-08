#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:53:42 2018

@author: guru

Script for Pei's interstitial scattering model for Cu in SnTe
"""

import numpy as np
from math import atan, pi
import matplotlib.pyplot as plt

eps = 83 
hbar = 1.054e-34
kb = 1.38e-23
MCu = 63.546
MSn = 118.71
MTe = 127.60
aSn = 6.3206
aCu = 6.3117
vs = 1640.075
vol = 31.256e-30

datafile = "/Users/guru/Documents/PDscattering/datafiles/CuSnTe_int.csv"
data = np.genfromtxt(datafile, dtype = float, delimiter = ',')
comp = data[:,0]
kappaL= data[:,1]
kappa0 = kappaL[0]
kappaf = kappaL[len(kappaL) -1]
GammaM = np.zeros(len(comp))
GammaInt = np.zeros(len(comp))
Gamma_a = np.zeros(len(comp))
Gamma = np.zeros(len(comp))
kL = np.zeros(len(comp))
u = np.zeros(len(comp))
kap_pure = np.zeros(len(comp))

#(mass) Scattering due to the Cu substituting on the Sn site
i=0
for x in comp:
    Mbar = MSn*(1-x) + MCu*x
    GammaM[i] = x*((MCu - Mbar)/Mbar)**2
    MbarInt = ((MSn + MTe)*(1-x) + (2*MCu + MTe)*x)/(2*(1-x) + 3*x)
    #MbarInt = x*MCu
    GammaInt[i] = x*((MCu)/MbarInt)**2
    a_bar = (aSn)*(1-x) + aCu*x
    Gamma_a[i] = eps*x*((aCu - a_bar)/a_bar)**2
    Gamma[i] = GammaM[i] + GammaInt[i] + Gamma_a[i]
    #Kappa expression
   # z = x/0.06
   # kap_pure[i] =  kappa0*(1-z) +  z*kappaf 
    prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)
    u[i] = (prefix*kappa0*Gamma[i])**(1/2)
    kL[i] = kappa0* atan(u[i])/u[i]
    i=i+1
#Plot
plt.scatter(comp, kappaL)
plt.plot(comp, kL)
