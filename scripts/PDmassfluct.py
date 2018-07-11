#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:27:08 2018
@author: guru

General script for mass fluctuation fitting
Journal paper: Thermal conductivity of gadolinium calcium silicate apatites
"""


"""
Read the csv file header to determine the number of components in the compound:  
"""
import numpy as np
from scipy.stats import binom 
import matplotlib.pyplot as plt
from math import pi as pi
"""
Constants
"""
kb = 1.38e-23
h = 1.054e-34


datafile = "/Users/guru/Documents/PDscattering/datafiles/AZnSb2_500K.csv"

#Parse header: length of header vector is number of site, and values are the site degenracies
with open(datafile, 'r') as f:
    stoich = f.readline().split("\t")
    stoich[-1] = stoich[-1][:-1]
    stoich = list(map(int, stoich))
    stoich = np.array(stoich)
    props = f.readline().split("\t")
    props[-1] = props[-1][:-1]
    mass = f.readline().split("\t")
    mass[-1] = mass[-1][:-1]
    mass = list(map(float, mass))
    mass = np.array(mass)
    subst = f.readline().split("\t")
    subst[-1] = subst[-1][:-1]
    subst = list(map(float, subst))
    subst = np.array(subst)
    nsites = len(stoich)
    natoms = sum(stoich)
    [debyeT, vol]  = list(map(float, props))
data = np.genfromtxt(datafile, dtype ='float', skip_header = 10, delimiter = ",")
comp = data[:,0]
kappa = (data[:,1])

vs = (kb/h)*debyeT*(vol/(6*pi**2))**(1/3)
"""
Plot the thermal conductivity data
"""
plt.scatter(comp, kappa)
plt.xlabel('Composition')
plt.ylabel('Thermal conductivity (W/mK)')

"""
Function: mass fluctuation parameter calculation with BFZ parameter
"""
def bfz(c,stoich, mass, subst):
    m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
    natoms = sum(stoich)
    gamma = 0
    for n in range(len(mass)):
        msite = subst[n]*c + mass[n]*(1-c)
        gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
        gamma = gamma + gamma_s
    return gamma

"""
Function: mass fluctuation parameter calculation wiht the PUC parameter 
"""

def puc(c,stoich, mass, subst):
    gammaPUC = 0
    m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
    #print(m_uc)
#    natoms = len(stoich)
    for n in range(len(mass)):
        gamma_s = 0
        for m in range(0,stoich[n]+1):
            #msite = subst[n]*c + mass[n]
            m_uc_nm = subst[n]*m + mass[n]*(stoich[n] - m) + np.dot(np.delete(stoich, n), np.delete(mass, n))*(1-c) + np.dot(np.delete(stoich, n), np.delete(subst, n))*(c)
           # print(m_uc_nm) 
            gamma_nm = binom.pmf(m, stoich[n], c)*(1-(m_uc_nm/m_uc))**2
            gamma_s = gamma_s + gamma_nm
        #print(gamma_s)
        gammaPUC = (gammaPUC+ gamma_s)
    return gammaPUC*natoms

"""
Function: mass fluctuation as mixture of end members due to electron counting
"""
def puc2(c,stoich, mass, subst):
    gamma = 0
    mH = np.dot(stoich, mass)
    mS = np.dot(stoich, subst)
    m_avg = mH*(1-c) + mS*c
    gamma = (c*(1-c)*((mH-mS)/m_avg)**2)
    return gamma*natoms

"""
Compute the thermal conductivity and plot
"""
kL = np.zeros(len(comp))
kL_PUC = np.zeros(len(comp))
kL_Yang = np.zeros(len(comp))
kL_PUC2 = np.zeros(len(comp))
prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)
i=0
crange = np.linspace(0.0001,0.99999,100)
gammaPUC = puc(0.9999, stoich, mass, subst)
for c in comp:
    gamma = bfz(c, stoich, mass, subst)
    gammaPUC = puc(c, stoich, mass, subst)
    gammaPUC2 = puc2(c,stoich, mass, subst)
    u = (prefix*gamma*kappa[0])**(1/2)
    uPUC = (prefix*gammaPUC*kappa[0])**(1/2)
    uPUC2 = (prefix*(gammaPUC2)*kappa[0])**(1/2)
    kL[i] =kappa[0]*np.arctan(u)/u
    kL_PUC[i] = kappa[0]*np.arctan(uPUC)/uPUC
    kL_PUC2[i] = kappa[0]*np.arctan(uPUC2)/uPUC2
    i=i+1
plt.plot(comp, kL, label = "mass fluctuation BFZ")
plt.plot(comp, kL_PUC, label = "mass fluctuation PUC")
plt.plot(comp, kL_PUC2, label = "mass fluctuation PUC limited")
plt.legend()

""" 
Yang's Model
"""
#kL_Yang = 1.2/((crange)*(1-crange))**(1/2) 
#plt.plot(crange, kL_Yang)
#plt.ylim([0,10])
#plt.savefig('Skutteruditepaper.pdf')