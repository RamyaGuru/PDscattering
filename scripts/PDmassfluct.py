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
import matplotlib as mpl


mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = '14'
"""
Constants
"""
kb = 1.38e-23
h = 1.054e-34


datafile = "/Users/guru/Documents/PDscattering/datafiles/CaZnMgSb_323K.csv"

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
    [debyeT, vol, grun, nu]  = list(map(float, props))
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
Function: mass fluctuation parameter calculation with the PUC parameter 
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
Function: includes vacancy scattering-- figure out how to generalize this
"""
def bfz_vac(c,stoich, mass, subst):
    m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
    natoms = sum(stoich)
    gamma = 0
    for n in range(len(mass)):
        if subst[n]==0 or mass[n]==0:
            msite = subst[n]*c + mass[n]*(1-c)
            gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*(-(abs(subst[n] - mass[n]))/msite - 2)**2
        else:
            msite = subst[n]*c + mass[n]*(1-c)
            gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
        gamma = gamma + gamma_s
    return gamma

"""
Function: vacancy scattering with fixed end members
"""
def puc_vac(c, stoich, mass, subst):
    gamma = 0
    mH = np.dot(stoich, mass)
    mS = np.dot(stoich, subst)
    m_avg = mH*(1-c) + mS*c
    gamma = (c*(1-c)*((-(abs(mS - mH))/m_avg - 2)**2))
    return gamma*natoms


"""
Compute the thermal conductivity and plot
"""
c_len = 100
kL = np.zeros(c_len)
kL_PUC = np.zeros(c_len)
kL_Yang = np.zeros(c_len)
kL_PUC2 = np.zeros(c_len)
kL_vac = np.zeros(c_len)
kap_pure = np.zeros(c_len)
prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/1640)
i=0
j=0
crange = np.linspace(comp[0],comp[len(comp)-1],100)
#gammaPUC = puc(0.9999, stoich, mass, subst)
for c in crange:
    gamma = bfz(c, stoich, mass, subst)
    gammaPUC = puc(c, stoich, mass, subst)
    gammaPUC2 = puc2(c,stoich, mass, subst)
    gammaVac = bfz_vac(c, stoich, mass, subst)
    kap_pure[j] = kappa[0]*(1-c) + kappa[len(kappa)-1]*(c)
    u = (prefix*gamma*kappa[0])**(1/2) 
    uPUC = (prefix*gammaPUC*kap_pure[j])**(1/2)
    uPUC2 = (prefix*(gammaPUC2)*kap_pure[j])**(1/2)
    uPUCvac = (prefix*(gammaVac)*kap_pure[j])**(1/2)
    if i==0:
        kL[i] =kappa[0]
        kL_PUC[i] = kappa[0]
        kL_PUC2[i] = kappa[0]
        kL_vac[i] = kappa[0]
    else:
        kL[i] =kappa[0]*np.arctan(u)/u
        kL_PUC[i] = kap_pure[j]*np.arctan(uPUC)/uPUC
        kL_PUC2[i] = kap_pure[j]*np.arctan(uPUC2)/uPUC2
        kL_vac[i] = kap_pure[j]*np.arctan(uPUCvac)/uPUCvac
    i=i+1
    j=j+1
plt.plot(crange, kL, label = "mass fluctuation BFZ")
plt.plot(crange, kL_PUC, label = "mass fluctuation PUC")
#plt.plot(crange, kL_PUC2, label = "mass fluctuation PUC limited")
#plt.plot(crange, kL_vac, label = "vacancy scattering")
plt.legend()

#%%
#"""
#Figure for the Skutterudite paper
#"""
#def skutFig():
#    pSkut = plt.figure()
#    ax = plt.axes(xlim = (0,1.0))
#    kL_Yang = (1/(kappa[0]*(1-crange) + kappa[len(kappa) - 1]*crange) + (1.2*((crange)*(1-crange))**(1/2)))**(-1)
#    plt.scatter(comp, kappa, color = 'xkcd:black', clip_on = False)
#    plt.plot(crange, kL_Yang, color = 'xkcd:green', label = 'Meisner et al.')
#    plt.plot(crange, kL_PUC, color = 'xkcd:dark purple', label = "This study")
#    plt.plot(crange, kL_PUC2, color = 'xkcd:dark purple', linestyle = "--", label = "U.C. Limited")
#    #plt.plot(crange, kL_vac, color = 'xkcd:black', linestyle = ":", label = 'mass + vacancy')
#    plt.fill_between(crange, kL_PUC, kL_PUC2, color= 'xkcd:pale purple', alpha = '0.5')
#    plt.ylabel('$\mathrm{\kappa_L (W/mK)}$')
#    plt.xlabel('$\mathrm{Composition (x)}$')
#    plt.ylim(0,10)
#    plt.xticks([0,0.2,0.4,0.6, 0.8, 1], ('0\n Co$_4$Sb$_{12}$', '0.2', '0.4', '0.6', '0.8', '1\n CeFe$_4$Sb$_{12}$'))
#    plt.legend(frameon =False)
#    #ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(6))
#    plt.savefig('Skutteruditepaper.png', dpi = 400, bbox_inches= 'tight')

#%%
"""
Figure for the Lanthanum paper
"""
#pLat = plt.figure()
#ax = plt.axes(xlim = (0,0.12))
#plt.scatter(comp, kappa, color = 'xkcd:black')
#plt.plot(crange, kL, color = 'xkcd:green', label = 'mass: Wang et al.')
#plt.plot(crange, kL_PUC, color = 'xkcd:purple', label = 'mass: This study')
#plt.plot(crange, kL_vac, color = 'xkcd:black', linestyle = '--', label = 'mass + bond')
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
##plt.xrange([0,0.12])
#ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(7))
#
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
#plt.xlabel(r'$\mathrm{Composition (x)}$')
#plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
#plt.legend(frameon = False)
#plt.savefig('Lan_paper.png', dpi = 400, bbox_inches= 'tight')


"""
Figure for the In2Te3-InSb Paper
"""
#pInTe3 = plt.figure()
#ax = plt.axes(xlim = (0,0.16))
#plt.scatter(comp, kappa, color = 'xkcd:black')
#plt.plot(crange, kL, color = 'xkcd:green', label = 'Pei et al.')
#plt.plot(crange, kL_PUC, color = 'xkcd:purple', label = 'mass: This study')
#plt.plot(crange, kL_vac, color = 'xkcd:black', linestyle = '--', label = 'mass + bond')
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
##plt.xrange([0,0.12])
#ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(5))
#
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
#plt.xlabel(r'$\mathrm{Composition (x)}$')
#plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
#plt.legend(frameon = False)
#plt.savefig('InTe_InSb.png', dpi = 400, bbox_inches= 'tight')

"""
Figure for the In2Te3-SnTe Paper
"""
#pInTe3 = plt.figure()
#ax = plt.axes(xlim = (0,0.12))
#plt.scatter(comp, kappa, color = 'xkcd:black')
#plt.plot(crange, kL, color = 'xkcd:green', label = 'Tan et al.')
#plt.plot(crange, kL_PUC, color = 'xkcd:purple', label = 'mass: This study')
#plt.plot(crange, kL_vac, color = 'xkcd:black', linestyle = '--', label = 'mass + bond')
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
##plt.xrange([0,0.12])
#ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(7))
#
##plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
#plt.xlabel(r'$\mathrm{Composition (x)}$')
#plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
#plt.legend(frameon = False)
#plt.savefig('InTe_SnTe.png', dpi = 400, bbox_inches= 'tight')

""" 
Yang's Model
"""
#kL_Yang = 1.2/((crange)*(1-crange))**(1/2) 
#plt.plot(crange, kL_Yang)
#plt.ylim([0,10])
#
#skutFig()
#plt.savefig('Skutteruditepaper_3curves.jpg')

#%%

"""
Plot for the CaZnSb Compounds
"""
pCZB = plt.figure()
ax = pCZB.add_axes([0,0,1,1])
ax = plt.axes(xlim = (0,1))
plt.scatter(comp, kappa, color = 'xkcd:black', clip_on=False)
plt.plot(crange, kL_PUC, color = 'xkcd:royal blue')
plt.plot(crange, kap_pure, color = 'xkcd:blue green')
#plt.plot(crange, kL_PUC, color = 'xkcd:purple', label = 'mass: This study')
#plt.plot(crange, kL_vac, color = 'xkcd:black', linestyle = '--', label = 'mass + bond')
#plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
#plt.xrange([0,0.12])
ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(6))

ax.text(0.5, 0.73, r'$\kappa_0$', verticalalignment = 'center', horizontalalignment = 'center', rotation =20, transform = ax.transAxes, color = 'xkcd:blue green', fontsize = 18)
ax.text(0.8, 0.3, r'$\kappa_0\frac{\mathrm{tan^{-1}}u}{u}$', verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes, color = 'xkcd:royal blue', fontsize = 18)

#plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
plt.xlabel(r'$\mathrm{(x)Ca(Mg_{1-x}Zn_x)_2Sb_2}$')
plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
plt.savefig('CZB_example.pdf', dpi = 400, bbox_inches= 'tight')