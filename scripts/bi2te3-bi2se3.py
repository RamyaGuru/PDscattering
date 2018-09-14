#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:44:38 2018

@author: guru

Site preference case for the Bi2Te3 and Bi2Se3
"""


from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt

directory  = "/Users/guru/Documents/PDscattering/datafiles/"

reg1Data= kLcurve(directory + 'bi2te3.csv')
reg2Data = kLcurve(directory + 'bi2sete2.csv')


reg1Fit = PropsM(directory + 'bi2te3.csv')
reg2Fit = PropsM(directory + 'bi2sete2.csv')

reg1Fit.data.comp = reg1Fit.data.comp*3
reg2Fit.data.comp = reg1Fit.data.comp*(3/2)


reg1Fit.KappaL()
reg2Fit.KappaL()

plt.scatter(reg1Data.comp, reg1Data.kappa)
plt.scatter(reg2Data.comp+ 0.3333333, reg2Data.kappa)

c1 = np.linspace(0.00000001, 0.3333333, 100)
c2 = np.linspace(0.3333333, 0.999999, 100)

plt.plot(c1, reg1Fit.kL)
plt.plot(c2 , reg2Fit.kL)

