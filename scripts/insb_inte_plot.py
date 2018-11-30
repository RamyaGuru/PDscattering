#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:19:44 2018

@author: guru

Fit the InSb and In2Te3 data without the "-2" term
"""

from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

directory = "/Users/guru/Documents/PDscattering/datafiles/"

insbData = kLcurve(directory + "InTe-InSb-300K.csv")
insbProps = PropsM(directory + "InTe-InSb-300K.csv")
insbProps.KappaL()

plt.scatter(insbData.comp, insbData.kappa)
crange = np.linspace(0.0000001,insbData.comp[len(insbData.comp)-1], 100)

plt.plot(crange, insbProps.kL)



"""
Figure for the In2Te3-InSb Paper
"""
pInTe3 = plt.figure()
ax = plt.axes(xlim = (0,0.16))
plt.scatter(insbData.comp, insbData.kappa, color = 'xkcd:black')
plt.plot(crange, kL, color = 'xkcd:dark purple', linestyle = '--', label = 'incorrect')
plt.plot(crange, kL_PUC, color = 'xkcd:purple', label = 'mass')
plt.plot(crange, kL_vac, color = 'xkcd:green', linestyle = '-.', label = 'mass + vacancy')
#plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
#plt.xrange([0,0.12])
ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(5))

#plt.xticks([0,0.02,0.04,0.06, 0.08, 0.1, 0.12], ('0.0', '0.02', '0.04', '0.06', '0.08', '0.1', '0.12'))
plt.xlabel(r'$\mathrm{Composition (x)}$')
plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
plt.legend(frameon = False)
plt.savefig('InTe_InSb.png', dpi = 200, bbox_inches= 'tight')