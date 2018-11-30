#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:35:28 2018

@author: guru

Phonon dispersion demos
"""

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = '14'
C = 1
m=1
M = 2
a =1
k = np.linspace(-pi/a, pi/a, 100)

#Acoustic branch

omega = np.sqrt(4*C/m)*np.abs(np.sin(k*a/2))

#Optical Branch

omegaO = np.sqrt(C*(m+M)/(m*M)+C*np.sqrt((M+m)**2/(M**2*m**2)- (4/(M*m))*np.square((np.sin(k*a/2)))))
omegaA = np.sqrt(C*(m+M)/(m*M)-C*np.sqrt((M+m)**2/(M**2*m**2)- (4/(M*m))*np.square((np.sin(k*a/2)))))

#Extended Acoustic Branch
a2 = 0.5
mbar = 1.5
k2 = np.linspace(-pi/a2, pi/a2, 100)
omegaA2 = np.sqrt(4*C/mbar)*np.abs(np.sin(k2*a2/2))
#omegaA2 = np.sqrt(C*(m+M)/(m*M)-C*np.sqrt((M+m)**2/(M**2*m**2)- (4/(M*m))*np.square((np.sin(k*a2/2)))))

frame1 = plt.gca()
ax = plt.axes(xlim = (-2*pi/a, 2*pi/a))
ax.yaxis.set_ticklabels([])
plt.xticks([-2*pi/a, 0, 2*pi/a], (r'$2\pi$/a', '0', r'$2\pi$/a'))
plt.plot(k2,omegaA2, color = 'xkcd:purple', linestyle = '--')
plt.plot(k, omegaA, color = 'xkcd:purple')
plt.plot(k,omegaO, color = 'xkcd:green')
plt.savefig('phonon_disp.jpg', dpi = 300, bbox_inches = 'tight')
plt.ylabel(r'$\omega(k)$')
plt.xlabel(r'k')