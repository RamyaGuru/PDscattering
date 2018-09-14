#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:34:42 2018

@author: guru

Test script for plotting normalized vacancy scattering all on one plot
"""
from PDscattering import kLcurve, PropsM
import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/guru/Documents/PDscattering/datafiles/"


#Get thermal conductivity data:
insbData = kLcurve(directory + "InTe-InSb-300K.csv")
snteData = kLcurve(directory + "InTe_SnTe.csv")
lanData = kLcurve(directory + "LaCoO_with_vac.csv")

#Get materials properties:
insbProps = PropsM(directory + "InTe-InSb-300K.csv")
snteProps = PropsM(directory + "InTe_SnTe.csv")
lanProps = PropsM(directory + "LaCoO_with_vac.csv")

#Normalize the thermal conductivity data:
insb_comp, insb_kappa = insbData.normalize()
snte_comp, snte_kappa = snteData.normalize()
lan_comp, lan_kappa = lanData.normalize()

#Compute the thermal conductivity 

insbProps.KappaL()
snteProps.KappaL()
lanProps.KappaL()

#Normalize the vacancy thermal conductivities

insb_kl = insbProps.normalize()
snte_kl = snteProps.normalize()
lan_kl = lanProps.normalize()

#Plot the vacancy curves

crange_insb = np.linspace(0.0000001,insbData.comp[len(insbData.comp)-1], 100)
crange_snte = np.linspace(0.0000001,snteData.comp[len(snteData.comp)-1], 100)
crange_lan = np.linspace(0.0000001,lanData.comp[len(lanData.comp)-1], 100)
plt.plot(crange_insb, insb_kl, color = 'xkcd:black')
plt.plot(crange_snte, snte_kl, color = 'xkcd:purple')
plt.plot(crange_lan, lan_kl, color = 'xkcd:green')

plt.scatter(insbData.comp, insb_kappa, label = 'InSb-In$_2$Te$_3$', color = 'xkcd:black')
plt.scatter(snteData.comp, snte_kappa, label = 'SnTe-In$_2$Te$_3$', color = 'xkcd:purple')
plt.scatter(lanData.comp, lan_kappa, label = 'LaCoO$_3$', color = 'xkcd:green')
plt.legend()


