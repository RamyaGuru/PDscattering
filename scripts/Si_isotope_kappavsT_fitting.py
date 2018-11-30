#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 18:01:53 2018

@author: guru

fit Gruneisen parameter from the kappa vs. T curves
using the Klemens/Callaway equations
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import pi as pi


#import data

TCond_Si= np.genfromtxt('/Users/guru/Documents/SnyderGroup/datafiles/TCond_Si.dat', skip_header=1)

comp=np.array([0,0.01, 0.05, 0.2])

#Rewrite as pandas dataframe

Tc_df= pd.DataFrame(data=TCond_Si, columns=['Temp', 'Si29_0', 'Si29_1', 'Si29_5', 'Si29_20', 'Si30_1','Si30_5', 'Si30_20'])

Tc_df.plot('Temp', ['Si29_0', 'Si29_1', 'Si29_5', 'Si29_20', 'Si30_1','Si30_5', 'Si30_20'], loglog = True)

Tc_df['Si30_0'] = Tc_df['Si29_0'] 
#Constants
k= 1.38e-23
Vol = 19.683e-30 #atomic volume
N=1/(Vol)
theta=498 #Debye Temperature
h= 1.054e-34
a= 3.867e-10 #Primitive unit cell parameters

#Plot the high temperature Umklapp scattering

Tc_df.plot('Temp', ['Si29_0', 'Si29_1', 'Si29_5', 'Si29_20', 'Si30_1','Si30_5', 'Si30_20'])


#Kappa vs. T equations-- High T Umklapp Scattering
M1 = 28
Trange = np.linspace(650,1000,45)
for M2 in [29,30]:    
    for c in comp:
        M = c*M2 + (1-c)*M1
        print(M)
        def kappaT(T, gamma): return (6*pi**2)**(2/3)/(4*pi**2)*(M*1.66e-27)/(T*Vol**(2/3)*gamma**2)*((theta*k/h)*(6*pi**2/Vol)**(-1/3))**3
        kCurve = 'Si' + str(M2) +'_'+ str(int(c*100))
        kappa = Tc_df.as_matrix(columns = [kCurve])
        temp = Tc_df.as_matrix(columns = ['Temp'])
        temp = temp.transpose()
        kappa = kappa.transpose()
        #print(np.shape(kappa))
        kappa = kappa[0,60:]
        temp = temp[0,60:]
        gamma, cov = curve_fit(kappaT, temp, kappa, bounds = (0, np.inf))
        print(gamma)
plt.figure()
plt.scatter(1/temp,kappa)
plt.xlim([0.001,0.002])
    
