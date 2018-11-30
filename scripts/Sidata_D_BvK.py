#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 13:44:04 2018
@author: guru
Test Debye model versus BvK 
"""

from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/guru/Documents/PDscattering/datafiles/"


##First, tried with the 800K data

Si_data = kLcurve(directory + "Si29_800K.csv")

Si_props = PropsM(directory + "Si29_800K.csv")

Si_props.KappaL()

crange_Si = np.linspace(0.0000001, Si_data.comp[len(Si_data.comp)-1], 100)

plt.plot(crange_Si, Si_props.kL, color = 'xkcd:purple', label = 'Debye model')
plt.scatter(Si_data.comp, Si_data.kappa, color = 'xkcd:dark purple')

Si_props.bvk_kappa()

plt.plot(crange_Si, Si_props.kL, color = 'xkcd:green', label = 'BvK model')
plt.scatter(Si_data.comp, Si_data.kappa, color = 'xkcd:dark purple')
plt.xlabel('Isotope Concentration %')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.legend()
plt.savefig('Si_disp_comp.jpg', dpi = 200, bbox_inches= 'tight')


##Now try with the 300K data 

Si_data = kLcurve(directory + "Si30_800K.csv")
Si_props = PropsM(directory + "Si30_800K.csv")

Si_props.KappaL()

crange_Si = np.linspace(0.0000001, Si_data.comp[len(Si_data.comp)-1], 100)

plt.figure()
plt.plot(crange_Si, Si_props.kL, color = 'xkcd:purple', label = 'Debye model')
plt.scatter(Si_data.comp, Si_data.kappa, color = 'xkcd:dark purple')

Si_props.bvk_kappa()

plt.plot(crange_Si, Si_props.kL, color = 'xkcd:green', label = 'BvK model')
plt.scatter(Si_data.comp, Si_data.kappa, color = 'xkcd:dark purple')
plt.xlabel('Isotope Concentration %')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.legend(frameon = False)

##Skutterudite Figure

Sk_data = kLcurve(directory + "FilledSkutterudite.csv")
Sk_props = PropsM(directory + "FilledSkutterudite.csv")

Sk_props.KappaL()

crange_Sk = np.linspace(0.0000001, Sk_data.comp[len(Sk_data.comp)-1], 100)

plt.figure()
plt.plot(crange_Sk, Sk_props.kL, color = 'xkcd:purple', label = 'Debye model')
plt.scatter(Sk_data.comp, Sk_data.kappa, color = 'xkcd:dark purple')
plt.ylabel('$\mathrm{\kappa_L (W/mK)}$')
plt.xlabel('$\mathrm{Composition (x)}$')
plt.ylim(0,10)
plt.xticks([0,0.2,0.4,0.6, 0.8, 1], ('0\n Co$_4$Sb$_{12}$', '0.2', '0.4', '0.6', '0.8', '1\n CeFe$_4$Sb$_{12}$'))

Sk_props.bvk_kappa()

plt.plot(crange_Sk, Sk_props.kL, color = 'xkcd:green', label = 'BvK model')
plt.scatter(Sk_data.comp, Sk_data.kappa, color = 'xkcd:dark purple')
plt.legend(frameon =False)
#ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(6))
plt.savefig('Skut_disp_comp.jpg', dpi = 200, bbox_inches= 'tight')