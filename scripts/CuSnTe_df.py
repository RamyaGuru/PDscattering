#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:26:23 2018

@author: guru
CuSnTe with the interstitials. Why isn;t the datafile misbehaving?
"""

import numpy as np
import matplotlib.pyplot as plt
from PDscattering import kLcurve, PropsM

directory = "/Users/guru/Documents/PDscattering/datafiles/"
data = kLcurve(directory + 'CuSnTe_df.csv')
props = PropsM(directory + 'CuSnTe_df.csv')

props.KappaL_bfz()
props.KappaL_puc()

crange = np.linspace(0.0000001, data.comp[len(data.comp)-1], 100)

plt.plot(crange, props.kLBFZ, color = 'xkcd:green')
plt.plot(crange, props.kLPUC, color = 'xkcd:purple')
plt.scatter(data.comp, data.kappa)

plt.figure()
plt.plot(crange, props.numbfz, color = 'xkcd:green')
plt.plot(crange, props.numpuc, color = 'xkcd:purple')