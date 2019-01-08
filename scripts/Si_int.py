#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 15:45:12 2018

@author: guru

Si interstitial test
"""

from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/guru/Documents/PDscattering/datafiles/"

#Get thermal conductivity data:
data = kLcurve(directory + "Si_int.csv")

#Calculate the mass difference parameter
props = PropsM(directory + "Si_int.csv")

#Si_props.KappaL_bfz()
props.KappaL_puc()
props.KappaL_bfz()

crange = np.linspace(0.0000001, data.comp[len(data.comp)-1], 100)

plt.plot(crange, props.kLPUC, color = 'xkcd:green')
plt.plot(crange, props.kLBFZ, color = 'xkcd:purple')

plt.scatter(data.comp, data.kappa)

#Compare to without the extra site:

data2 =kLcurve(directory + "Si_int2.csv")

props2 = PropsM(directory + "Si_int2.csv")

props2.KappaL_puc()
props2.KappaL_bfz()
plt.plot(crange, props2.kLPUC, linestyle = ':', color = 'xkcd:green')
plt.plot(crange, props2.kLBFZ, linestyle = ':', color = 'xkcd:purple')
plt.scatter(data2.comp, data2.kappa)

