#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:17:50 2018

@author: guru
"""
"""
Modelling of the thermal conducitvity versus temperature for silicon and its isotopes. Uses the Goldsmid formulation. 

Need to introduce factor "n" represeting the number of atoms in the unit cell which produces the best fit
"""
from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mtplt

mtplt.rcParams['font.size'] = '16'
mtplt.rcParams['font.family'] = 'sans-serif'
mtplt.rcParams['lines.linewidth'] = 2
mtplt.rcParams['font.weight'] = 'bold'
mtplt.rcParams['axes.xmargin'] = 0.1
mtplt.rcParams['axes.titlesize'] = 'Large'
mtplt.rcParams['axes.formatter.useoffset'] = False


TCond_Si= np.genfromtxt('/Users/guru/Documents/SnyderGroup/datafiles/TCond_Si.dat', skip_header=1)

comp=np.array([0,0.01, 0.05, 0.2])

#Rewrite as pandas dataframe

Tc_df= pd.DataFrame(data=TCond_Si, columns=['Temp', 'Pure', 'Si29_1', 'Si29_5', 'Si29_20', 'Si30_1','Si30_5', 'Si30_20'])

#Generate arrays for Si29 and Si30 at 800K as this is above the Debye temperature for Silicon (640K)
temp= 800#choose a temperature above the Debye temperature to perform the anaylsis

row = Tc_df.loc[Tc_df['Temp']==temp]
Tc_Si30= row.as_matrix(columns=['Pure' ,'Si30_1', 'Si30_5', 'Si30_20'])
Tc_Si30=Tc_Si30.transpose()
Tc_Si29= row.as_matrix(columns=['Pure' ,'Si29_1', 'Si29_5', 'Si29_20'])
Tc_Si29=Tc_Si29.transpose()
p30= plt.scatter(comp, Tc_Si30)
p29= plt.scatter(comp, Tc_Si29)

#Fit with the Goldsmid formulation

# =============================================================================
# Silicon parameters
# Debye temperature: 640 K
# Volume per atom: 2.002e-29 m^3
# atoms per unit cell: 8  
# =============================================================================

#Constants
k= 1.38e-23
V= 2.0012875875e-29#atomic volume
#n= 8 #atoms per cubic unit cell
#N=4*6.2459788778e27 #Concentration of primitive unit cells (1/(number of atoms in primitive unitcell)*(atomic volume))
N=1/(V)
theta=531.73
h= 1.054e-34
a= 5.43e-10




m0=28
#mi=29
#mi2=30

k0_Si30= Tc_Si30[0];
k0_Si29= Tc_Si29[0];

#c is the composition, plot between 0 and 0.25

c= np.linspace(0.001,0.25, 100)
#c=0.001
index=1;
size=1000
mass= np.array([29,30])
eps_Si=np.zeros((100,))
#n=0
z=6
for mi in mass:
    if index ==1:
    #Atom scattering parameters
        eps_Si_n=c*(1-(mi/(mi*c + m0*(1-c))))**2 + (1-c)*(1-(m0/(mi*c + m0*(1-c))))**2
    elif index ==2:
    #Abridged formula from Abeles 
        mavg= mi*c + m0*(1-c)
        eps_Si_n= c*(1-c)*((mi-m0)/mavg)**2
    elif index ==3:
    #Unit cell scattering parameter 1 (Riley)
        eps_Si_n=2*c*(1-((mi+m0)/(2*(mi*c + m0*(1-c)))))**2 + (1-2*c)*(1-((m0+m0)/(2*(mi*c + m0*(1-c)))))**2
    elif index ==4:
    #Unit cell scaterring parameter 2 (Riley)
        eps_Si_n=c*(1-((mi+mi)/(2*(mi*c + m0*(1-c)))))**2 + (1-c)*(1-((m0+m0)/(2*(mi*c + m0*(1-c)))))**2
    elif index ==5:
    #2 atom unit cell with Klemens model (from Yinglu paper)
        cii= c**2;
        ci0=c*(1-c)*2;
        c00=(1-c)**2;
        mavg= (mi+mi)*cii + (mi+m0)*ci0 + (m0+m0)*c00
        mavg= 2*(mi*c + m0*(1-c))
        eps_Si_n=cii*(1-((mi+mi)/mavg))**2 + ci0*(1-((mi+m0)/mavg))**2 + c00*(1-((m0+m0)/mavg))**2
    elif index ==6:
    #Unit cell scattering parameter 4 (Ramya)
        eps_Si_n=c*(1-((mi+mi)/((mi+mi)*c + (m0+m0)*(1-c))))**2 + (1-c)*(1-((m0+m0)/((mi+mi)*c + (m0+m0)*(1-c))))**2
    elif index ==7:
    #Large unit cell containing four atoms
        eps_Si_n=c*z*(1-((mi+m0*(z-1))/((mi+m0*(z-1))*c*z + (m0*z)*(1-c*z))))**2 + (1-c*z)*(1-((m0*z)/((mi+m0*(z-1))*c*z + (m0*z)*(1-c*z))))**2 
    elif index ==8:
        #4 atom unit cell with appropriate probabilities
        ciiii=c**4;
        ciii=c**3*(1-c)*4
        cii=c**2*(1-c)**2*6
        ci= c*(1-c)**3*4
        c0=(1-c)**4
        mavg= mi*c + m0*(1-c)
        eps_Si_n=ciiii*(1-((mi)/mavg))**2 + ciii*(1-(((mi*3+m0)/4)/mavg))**2 + cii*(1-(((mi*2+m0*2)/4)/mavg))**2 + ci*(1-(((mi+m0*3)/4)/mavg))**2 + c0*(1-(m0/mavg))**2
    elif index==9:
        #virtual crystal approximation: average over the mass of the primitive unit cells 
        cii= c**2;
        ci0=c*(1-c)*2;
        c00=(1-c)**2;
        mavg = mi*c + m0*(1-c)
        eps_Si_n=cii*(1-((mi)/mavg))**2 + ci0*(1-(((mi+m0)/2)/mavg))**2 + c00*(1-((m0)/mavg))**2
 
        
    eps_Si = np.vstack((eps_Si, eps_Si_n))  

        
    
    
omega_Si=np.zeros((100,))
#vals= [row[1] for row in eps_Si]
for r in [1,2]:
    omega_Si_n= ((pi*h/(2*k**2))*(6*pi**2/V)**(2/3)*k0_Si29*(eps_Si[r,:])/(N*theta))**(1/2)
    omega_Si= np.vstack((omega_Si, omega_Si_n))
    print(omega_Si.shape)
#omega_Si30= ((pi*h/(2*k**2))*(6*pi**2/V)**(2/3)*k0_Si30*(eps_Si30)/(N*theta))**(1/2)
#omega_Si29= ((pi*h/(2*k**2))*(6*pi**2/V)**(2/3)*k0_Si29*(c*((29-(29*c+ 28*(1-c)))**2/(29*c+ 28*(1-c))**2)+ (1-c)*((28-(29*c+ 28*(1-c)))**2/(29*c+ 28*(1-c))**2))/(N*theta))**(1/2)

kl_Si29= k0_Si29 * (np.arctan(omega_Si[1,:])/omega_Si[1,:])
kl_Si30=k0_Si30 * (np.arctan(omega_Si[2,:])/omega_Si[2,:])
#print(omega_Si[1,:])
#print(kl_Si29)
font = {'family' : 'normal',
        'size'   : 16}

mtplt.rc('font', **font)
mtplt.rcParams['lines.linewidth'] = 3
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
p1 = plt.scatter(comp*100, Tc_Si29, color= 'xkcd:dark purple', label='Si$_{29}$ data')
p2 = plt.scatter(comp*100, Tc_Si30, color ='xkcd:dark green', label='Si$_{30}$ data')
p3, = plt.plot(c*100, kl_Si29, color='xkcd:purple', label='Si$_{29}$ Klemens')
p4, = plt.plot(c*100, kl_Si30, color='xkcd:green', label='Si$_{30}$ Klemens')

#plt.legend()
#plt.legend([p1,p2, p3, p4,],['Si29 data', 'Si30 data', 'Si29_fit', 'Si30_fit'])
#plt.xlim([0,1])
#plt.gcf().subplots_adjust(bottom=0.15)
#plt.gcf().subplots_adjust(left=0.15)
#plt.xlim([0,0.1])
#plt.ylim([55,58])


ax.text(0.5,0.83,'Si-29', horizontalalignment = 'center', verticalalignment = 'center', rotation =-10, transform = ax.transAxes, color = 'xkcd:dark purple', fontsize = 16)
ax.text(0.5, 0.42, 'Si-30', horizontalalignment = 'center', verticalalignment = 'center', rotation = -36, transform = ax.transAxes, color = 'xkcd:dark green', fontsize = 16)

plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'$\mathrm{Si^{\prime}_{x}Si_{1-x}}$')
plt.savefig('isotope_scattering_site_800K.pdf', dpi=300, bbox_inches = 'tight')