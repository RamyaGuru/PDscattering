#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 11:33:28 2018

@author: guru


Skutterudite Test file
"""

import numpy as np
import matplotlib.pyplot as plt
from PDscattering import kLcurve, PropsM
import os 

directory = "/Users/guru/Documents/PDscattering/SkuttFilling/datafiles/"

for file in os.listdir(directory):
    if file.endswith(".dat"):
        data = kLcurve(directory + file)
        props = PropsM(directory + file)
        props.KappaL()
        # Plot curves
        crange = np.linspace(0.0000001,data.comp[len(data.comp)-1], 100)
        plt.figure()
        plt.plot(crange, props.kL)
        plt.scatter(data.comp, data.kappa)
        plt.xlabel('Filling Fraction (x)')
        plt.ylabel(r'$\kappa_L$ (W/m/K)')
        plt.savefig(file + '.jpg', dpi = 300, bbox_inches = 'tight')
        np.savetxt(file + '.txt', props.kL, delimiter = ",")

