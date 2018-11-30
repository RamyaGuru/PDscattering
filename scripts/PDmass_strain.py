#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:23:36 2018

@author: guru

PD mass and strain fluctuation scattering



Note to self: need to adjust for the thermal conductivity of the pure materials

how to separately treat the U-curve?
"""

from PDmassfluct import bfz, puc, skutFig
import numpy as np
import matplotlib.pyplot as plt
from math import pi as pi 
from scipy.stats import binom
from scipy.optimize import curve_fit
import matplotlib as mpl
#from lmfit import minimize, Parameters

mpl.rcParams['font.size'] = '14'
#%%
"""
Constants
"""
kb = 1.38e-23
#debyeT = 386.554
#vol = 12.6e-30
h = 1.054e-34
eps = 3.5
eps2= 0.21
datafile = "/Users/guru/Documents/PDscattering/datafiles/bi2sete2.csv"

#%%
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
    rad = f.readline().split("\t")
    rad[-1] = rad[-1][:-1]
    rad = list(map(float, rad))
    rad = np.array(rad)
    radsub = f.readline().split("\t")
    radsub[-1] = radsub[-1][:-1]
    radsub = list(map(float, radsub))
    radsub = np.array(radsub)
    nsites = len(stoich)
    natoms = sum(stoich)
    lat = f.readline().split("\t")
    lat[-1] = lat[-1][:-1]
    lat = list(map(float, lat))
    lat = np.array(lat)
    latsub = f.readline().split("\t")
    latsub[-1] = latsub[-1][:-1]
    latsub = list(map(float, latsub))
    latsub = np.array(latsub)
    G = f.readline().split("\t")
    G[-1] = G[-1][:-1]
    G = list(map(float, G))
    G = np.array(G)
    Gsub = f.readline().split("\t")
    Gsub[-1] = Gsub[-1][:-1]
    Gsub = list(map(float, Gsub))
    Gsub = np.array(Gsub)
    [debyeT, vol, grun, nu]  = list(map(float, props))
data = np.genfromtxt(datafile, dtype ='float', skip_header = 10, delimiter = ",")
comp = data[:,0]
kappa = data[:,1]
vs = (kb/h)*debyeT*(vol/(6*pi**2))**(1/3)
c_len = 100

kL_BFZ = np.zeros(c_len)
kL_PUC = np.zeros(c_len)
kL_MV = np.zeros(c_len)
kL_MG = np.zeros(c_len)
kL_MV2 = np.zeros(c_len)
kL_MVPUC = np.zeros(c_len)
kL_G = np.zeros(c_len)
kL_VPUC2 = np.zeros(c_len)

print(radsub)
#%%
"""
Functions for the strain Scattering Parameter
"""
#Yang function uses exact inputs from the Yang model
def yangmod(c, stoich, rad, radsub):
    natoms = sum(stoich)
    r_uc = (np.dot(stoich, radsub)*c + np.dot(stoich, rad)*(1-c))/natoms
    natoms = sum(stoich)
    gamma = 0
    for n in range(len(mass)):
        rsite = radsub[n]*c + rad[n]*(1-c)
        gamma_s = (stoich[n]/natoms)*c*(1-c)*(rsite/(r_uc))**2*((rad[n] - radsub[n])/rsite)**2
        gamma = gamma + gamma_s
    return gamma

#YangMod uses the radii instead while making the monatomic lattice approximation
def yang(c,stoich, rad, radsub):
    natoms = sum(stoich)
    m_uc = (np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c))/natoms
    natoms = sum(stoich)
    gamma = 0
    for n in range(len(mass)):
        rsite = radsub[n]*c + rad[n]*(1-c)
        msite = subst[n]*c + mass[n]*(1-c)
        gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc))**2*((rad[n] - radsub[n])/rsite)**2
        gamma = gamma + gamma_s
    return gamma


#PUC: Stat mech interpretation using the change in the lattice parameter
    #of the primitive unit cell to describe the system.
    #Apply Vegard's law for full statistical mechanical handling
def pucV(c,stoich, lat, latsub):
    l_uc = latsub*c + lat*(1-c)
    gamma = c*(1-c)*((lat-latsub)/l_uc)**2            
    return gamma

def pucV2(c,stoich, lat, latsub):
    l_uc = latsub*c + lat*(1-c)
    for n in range(len(mass)):
        for m in range(0,stoich[n]+1):
            lat_eff = lat*(1-(m/stoich[n])) + latsub*(m/stoich[n])
            gamma = binom.pmf(m,stoich[n],c)*(1-binom.pmf(m,stoich[n],c))*(1-((lat-lat_eff)/l_uc))**2          
    return gamma

"""
Functions for force constant fluctuation parameters 
"""

#Solid solution between two types of unit cells
def pucG(c,stoich, G, Gsub):
    G_uc = G
    gamma = c*(1-c)*(1-((G-Gsub)/G_uc))**2
    return gamma

#Vegard's law for stat. mech. model:


"""
Function to plot the thermal conductivity data
"""
def pltdata():
    plt.scatter(comp, kappa)
    plt.xlabel('Composition')
    plt.ylabel('Thermal conductivity (W/mK)')

#%%
"""
Function for computing kappa from gamma
"""
def TCond(gamma, x):
    kap_pure = kappa[0]*(1-x) + kappa[len(kappa) - 1]*x
    u = (prefix* gamma* kap_pure)**(1/2)
    kL = kap_pure *(np.arctan(u)/u)
    return kL
i=0
prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)

"""
Functions for fitting epsilon
"""
def fiteps(Mfunc, Sfunc, inputM, inputS, comp, kappa):
    def gammatot(c,epsilon): 
        gamma = Mfunc(c, *inputM) + epsilon*Sfunc(c, *inputS) 
        kL = TCond(gamma, c)
        return kL
    eps, cov = curve_fit(gammatot, comp, kappa, bounds = (0,np.inf))
    return eps
##%%
#"""
#Computing epsilon
#"""
#Gfluct = abs(G-Gsub)/Gsub
#latfluct = abs(lat - latsub)/latsub
#C = Gfluct/latfluct
#rat = Gsub/G
#grun = 0.95 #3.13 # gruneisen for SnSe from Xiao2016 
#nu = 0.193 #0.21806 #Voigt Posson ratio for SnSe from MatProj
#conv= (1+nu)*rat/(2+rat+nu*(rat-4))
epscalc = 12*(conv*((1/6)**(1/2)*C - 4.2*grun*(2/3)**(1/2)))**2
#epscalc = 2*((C-6.4*grun)*((1+nu)*rat/(2+rat+nu*(rat-4))))**2
eps_noG = 2*(6.4*grun*((1+nu)/(3*(1-nu))))**2
#print('Calculated value for epsilon is: %4.4f' % epscalc)
print('Calculated value of epsilon (no G fluctuation): %4.4f' % eps_noG)
""" 
Function call to get the value of epsilon 
"""
inputM = (stoich, mass, subst)
inputS = (stoich, rad, radsub)
pltdata()
epsMVyang = fiteps(bfz, yang, (stoich, mass, subst), (stoich, rad, radsub), comp, kappa)
print('Epsilon for BFZ mass and Yang volume: %4.4f' % epsMVyang)
epsMVyangmod = fiteps(bfz, yangmod, (stoich, mass, subst), (stoich, rad, radsub), comp, kappa)
print('Epsilon for BFZ mass and Yangmod: %4.4f' % epsMVyangmod)
#epsMVpuc= fiteps(bfz, pucV, (stoich, mass, subst), (stoich, lat, latsub), comp, kappa)
#print('Epsilon for BFZ mass and PUC vol: %4.4f' % epsMVpuc)
for c in np.linspace(0.0001,0.9999,100):
    inputM = (stoich, mass, subst)
    inputS = (stoich, rad, radsub)
    gammaM_bfz = bfz(c,stoich, mass, subst)
    gammaM_puc = puc(c,stoich, mass, subst)
    gammaV_yang = yang(c,stoich, rad, radsub)
    gammaV_yangmod = yangmod(c, stoich, rad, radsub)
    #gammaV_puc = pucV(c, stoich, lat, latsub)
    #gammaV_puc2 = pucV2(c, stoich, lat, latsub)
    #gammaG = pucG(c, stoich, G, Gsub)
    gammaMVyang = gammaM_bfz + epsMVyang* gammaV_yang
    gammaMVyangmod = gammaM_bfz + epsMVyangmod*gammaV_yangmod
    #gammaMVpuc = gammaM_bfz + epsMVpuc*gammaV_puc
    #gammaMG = gammaM_bfz + gammaG
    #gamma1 = gammaM_bfz + gammaV_yang
    #gamma2 = gammaM_bfz + gammaV_yangmod
    #gamma3 = gammaM_bfz + gammaV_puc
#    uBFZ = (prefix*gammaM_bfz*kappa[0])**(1/2)
#    uPUC = (2*prefix*gammaM_puc*kappa[0])**(1/2)
#    uV_PUC = (prefix*gammaV_puc*kappa[0])**(1/2)
#    uV_PUC2 = (prefix*gammaV_puc2*kappa[0])**(1/2)
#    uG = (0.6*prefix*gammaG*kappa[0])**(1/2)
#    u1 = (eps*prefix*gamma1*kappa[0])**(1/2)
    #u2 = (eps*prefix*gamma2*kappa[0])**(1/2)
    #u3 = (2*eps2*prefix*gamma3*kappa[0])**(1/2)
#    kL_BFZ[i] =kappa[0]*np.arctan(uBFZ)/uBFZ
#    kL_PUC[i] = kappa[0]*np.arctan(uPUC)/uPUC
#    kL_VPUC2[i] = kappa[0]*np.arctan(uV_PUC2)/uV_PUC2
#    kL_1[i] = kappa[0]*np.arctan(u1)/u1
    #kL_2[i] = kappa[0]*np.arctan(u2)/u2
    #kL_3[i] = kappa[0]*np.arctan(u3)/u3
#    kL_G[i] = kappa[0]*np.arctan(uG)/uG
    kL_MV[i] = TCond(gammaMVyang, c)
    kL_MV2[i] = TCond(gammaMVyangmod, c)
    #kL_MG[i] = TCond(gammaMG, c)
    #kL_G[i] = TCond(gammaG, c)
    #kL_MVPUC[i] = TCond(gammaMVpuc, c)
    i = i+1  
    
#%%
crange = np.linspace(0.0001,0.9999,100)

"""
Plot mass fluctuation alone
"""
#plt.figure()
#pltdata()
#plt.plot(comp, kL_BFZ, label = 'mass BFZ')
#plt.plot(comp, kL_PUC, label = 'mass PUC')

"""
Plot of volume flucatuation alone

"""
plt.figure()
pltdata()
plt.plot(crange, kL_MV)

"""
Plot mass and volume fluctuation alone
"""
#plt.figure()
#pltdata()
#plt.plot(comp, kL_MV, label = 'mass vol Yang')
#plt.plot(comp, kL_MG, label = 'mass G')
#plt.plot(comp, kL_MV2, label = 'mass vol Yangmod')
#plt.plot(comp, kL_MVPUC, label  = 'mass vol PUC')
#plt.plot(comp, kL_G, label = 'G PUC')
##plt.plot(comp, kL_)
#plt.legend()

#%%
"""
Plot of strain scattering for paper
"""

skutFig()
plt.plot(crange, kL_MVPUC, label  = 'mass vol PUC', color = 'xkcd:black', linestyle = ':')
plt.savefig('FilledSkutterudite_v2.png', dpi = 200, bbox_inches = 'tight')

#%%

#plt.plot()



#def yangmod(stoich, rad, radsub, c):


#def pucVol(stoich, rad, radsub, c):
    