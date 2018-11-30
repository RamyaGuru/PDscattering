#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 13:37:06 2018

@author: guru

plotting group velocities BvK vs. Debye
"""
import numpy as np
from math import pi as pi
import matplotlib.pyplot as plt
import scipy.integrate as integrate

vol = 19.68e-30
vs = (kb/h)*531.73*(vol/(6*pi**2))**(1/3)
kmax = (6*pi**2/vol)**(1/3) 
k = np.linspace(0, kmax, 100)
kap_pure = 57.78
kb = 1.38e-23
gamma = 0.5*(1/28.5)**2
#BvK
vg3BvK = (vs*(np.cos((pi/2)* (k/kmax))))**3

#Debye
vg3D = vs**3+0*k

plt.plot(k, vg3D, color = 'xkcd:green', label= 'Debye')
plt.plot(k,vg3BvK, color = 'xkcd:purple', label = 'BvK')

plt.xlabel(r'$k$')
plt.ylabel(r'$v_g^3$')
plt.legend()
plt.savefig('vg3_comp.png', dpi = 300, bbox_inches = 'tight')

plt.figure()

#Arctan plots: strong PD scattering
omegaD = kmax*vs
omega = np.linspace(0,omegaD)

def vg3Deb(k):
    vg3 = vs**3+0*k
    return vg3
def vg3BvK(k):
    vg3 = (vs*(np.cos((pi/2)* (k/kmax))))**3
    return vg3
vg3_avgD, other = integrate.quad(vg3Deb, 0, kmax)
vg3_avgBvK, other = integrate.quad(vg3BvK, 0, kmax)

#BvK
a = vol/(4*pi)*gamma
bBvK = kb/(2*pi**2*kap_pure)*vg3_avgBvK
bD = kb/(2*pi**2*kap_pure)*vg3_avgD

#BvK
atanIntBvK = 1/(1+ a*omega**2/bBvK)
atanIntD = 1/(1+ a*omega**2/bD)
plt.plot(omega, atanIntBvK, color = 'xkcd:purple', label = 'BvK')

plt.xlabel(r'$\omega$')
plt.ylabel(r'$1/(1+(a\omega^2/b)$')
plt.savefig('integrandForm.png', dpi = 300, bbox_inches = 'tight')
#plt.plot(omega, atanIntD, color = 'xkcd:green', label= 'Debye')