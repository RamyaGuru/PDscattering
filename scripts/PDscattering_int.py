#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:59:08 2018

@author: guru

Class defiitons:
classes:
    kLcurve: attributes will include kappa and comp
            methods will include init and normalize
    matPropMass: attributes will include mass, subst...
                 methods will include calculation of gamma.. for special cases too
    matPropMass+Struct: attributes will include mass, subst, vol, grun...
                        methods will include calculation of gamma and epsilons
    Gamma: attributes will include 
"""

import numpy as np
from math import pi as pi
from scipy.stats import binom
from math import sin, cos, atan
import scipy.integrate as integrate

import matplotlib.pyplot as plt
"""
Constants
"""
kb = 1.38e-23
h = 1.054e-34


class kLcurve:
    '''
    Args: filepath (str)
    
    Attributes:
        .kappa -> thermal conductivity data
        .comp -> site concnetration of defects
    Methods:
        .normalize(kappa) -> normalize the data so that thermal conductivity varies between 0 and 1
        .plotdata(kappa, c) -> nicely formatted graph for scatter plot of data
    '''

    def __init__(self, fileID):
        data = np.genfromtxt(fileID, dtype ='float', skip_header = 10, delimiter = ",")
        cols = data.shape[1]
        self.comp = data[:,0]
        self.kappa = (data[:, cols-1])
    
    def normalize(self):
        kappa = self.kappa
        comp = self.comp
        return comp/(max(comp)), kappa/(max(kappa))
    
    
    
    
class PropsM:
    
    def __init__(self,fileID):   
        with open(fileID, 'r') as f:
            stoich = f.readline().split("\t")
            stoich[-1] = stoich[-1][:-1]
            stoich = list(map(int, stoich))
            self.stoich = np.array(stoich)
            props = f.readline().split("\t")
            props[-1] = props[-1][:-1]
            self.props = props
            mass = f.readline().split("\t")
            mass[-1] = mass[-1][:-1]
            mass = list(map(float, mass))
            self.mass = np.array(mass)
            subst = f.readline().split("\t")
            subst[-1] = subst[-1][:-1]
            subst = list(map(float, subst))
            self.subst = np.array(subst)
#            for n in range(len(subst)):
#                strlist = subst[n].split(",")
#                subst[n] =[float(i) for i in strlist]
            self.subst = subst
            self.nsites = len(stoich)
            self.natoms = sum(stoich)
            [self.debyeT, self.vol] = [list(map(float, props))[i] for i in [0,1]]
            self.vs = (kb/h)*self.debyeT*(self.vol/(6*pi**2))**(1/3)
            self.data = kLcurve(fileID)
            if self.data.shape[1] != len(stoich):
                raise ValueError('give composition for each impurity')
    
    """
    Function: mass fluctuation parameter calculation with BFZ parameter
    """
    def bfz(self, c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
        natoms = sum(stoich)
        gamma = 0
        for n in range(len(mass)):
            msite = subst[n]*c + mass[n]*(1-c)
            #gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
            gamma_s = (stoich[n]/natoms)*(msite/(m_uc/natoms))**2*(c*((subst[n] - msite)/msite)**2+(1-c)*((mass[n] - msite)/msite)**2)
            print(msite)
            gamma = gamma + gamma_s
            #print("For site %d, the gamma is: %4.8f" % (n, gamma_s))
        self.gamma = gamma
            
            

 '''
 Function: PUC 
 '''


    def puc(self,c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
        natoms = sum(stoich)
        sblatt = len(stoich)
        #full_p = 0
        #a = 0
        pUC = [1]
        pUC_pregen = [1]*sblatt
        m_def_pregen = [0]*sblatt
        m_def = [0]
        GammaTot = [0]
        func = [0]*sblatt
        #Have a counter for which function call I am in
        def recloop(a):
            func[a] +=1
            for m in range(stoich[a]+1):
                pUC[0] = pUC_pregen[a] #should subtract based on the function call number... 
                m_def[0] = m_def_pregen[a] 
                pUC[0] = pUC[0]*binom.pmf(m,stoich[a],c)
                m_def[0] = m_def[0] + subst[a]*m +mass[a]*(stoich[a] - m)
                if a+1<sblatt:
                    pUC_pregen[a+1] = pUC[0]
                    m_def_pregen[a+1] = m_def[0]
                    recloop(a+1)
                else:
                    gamma_nm = (1-(m_def[0])/m_uc)**2
                    GammaTot[0] = GammaTot[0] + pUC[0]*gamma_nm
        recloop(0)
        self.gamma = GammaTot[0]*natoms
        
    def puc2(self,c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
        natoms = sum(stoich)
        sblatt = len(stoich)
        #full_p = 0
        #a = 0
        pUC = [1]
        pUC_pregen = [1]*sblatt
        m_def_pregen = [0]*sblatt
        m_def = [0]
        GammaTot = [0]
        func = [0]*sblatt
        #Have a counter for which function call I am in
        def recloop(a):
            func[a] +=1
            for m in range(stoich[a]+1):
                pUC[0] = pUC_pregen[a] #should subtract based on the function call number... 
                m_def[0] = m_def_pregen[a] 
                pUC[0] = pUC[0]*binom.pmf(m,stoich[a],c)
                m_def[0] = m_def[0] + subst[a]*m +mass[a]*(stoich[a] - m)
                if a+1<sblatt:
                    pUC_pregen[a+1] = pUC[0]
                    m_def_pregen[a+1] = m_def[0]
                    recloop(a+1)
                else:
                    gamma_nm = (1-(m_def[0])/m_uc)**2
                    GammaTot[0] = GammaTot[0] + pUC[0]*gamma_nm
        recloop(0)
        self.gamma = GammaTot[0]*natoms        
                        
          

        
    """
    Function: mass fluctuation as mixture of end members due to electron counting
    """
    def puc2(self, c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        natoms = sum(stoich)
        gamma = 0
        mH = np.dot(stoich, mass)
        mS = np.dot(stoich, subst)
        m_avg = mH*(1-c) + mS*c
        gamma = (c*(1-c)*((mH-mS)/m_avg)**2)
        self.gamma = gamma*natoms
        
    
    """
    Function: includes vacancy scattering-- figure out how to generalize this
    """
    def bfz_vac(self, c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        natoms = sum(stoich)
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
        gamma = 0
        for n in range(len(mass)):
            if subst[n]==0:
                msite = subst[n]*c + mass[n]*(1-c)
                gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*(-(3*(abs(subst[n] - mass[n]))/msite))**2
            else:
                msite = subst[n]*c + mass[n]*(1-c)
                gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
            gamma = gamma + gamma_s
        self.gamma = gamma


    """
    Function: calculate the thermal conductivity from gamma-- BFZ method
    """
    def KappaL_bfz(self): #for u-curves... considers both end-members
        num = sum(self.stoich)
        comp = self.data.comp
        vol = self.vol
        vs = self.vs
        prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)
        c_len = 100
        kL = np.zeros(c_len)
        crange = np.linspace(comp[0], comp[len(comp)-1], 100)
        i = 0
        for c in crange:
            self.bfz(c)
            #self.bfz_vac(c)
            if comp[len(comp)-1]>0.9:
                kap_pure = self.data.kappa[0]*(1-c) + self.data.kappa[len(self.data.kappa)-1]*c
            else:
                kap_pure = self.data.kappa[0]
            u = (prefix*self.gamma*kap_pure)**(1/2)
            if i==0:
                kL[i] = self.data.kappa[0]
            else:
                kL[i] = kap_pure*np.arctan(u)/u
            i = i+1
        self.kLBFZ = kL
    
    """
    Function: calculate the thermal conductivity from gamma-- PUC method
    """
    def KappaL_puc(self): #for u-curves... considers both end-members
        num = sum(self.stoich)
        comp = self.data.comp
        vol = self.vol
        vs = self.vs
        prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)
        c_len = 100
        kL = np.zeros(c_len)
        crange = np.linspace(comp[0], comp[len(comp)-1], 100)
        i = 0
        for c in crange:
            self.puc(c)
            #self.bfz_vac(c)
            if comp[len(comp)-1]>0.9:
                kap_pure = self.data.kappa[0]*(1-c) + self.data.kappa[len(self.data.kappa)-1]*c
            else:
                kap_pure = self.data.kappa[0]
            u = (prefix*self.gamma*kap_pure)**(1/2)
            print(self.gamma)
            if i==0:
                kL[i] = self.data.kappa[0]
            else:
                kL[i] = kap_pure*np.arctan(u)/u
            print(u)
            i = i+1
        self.kLPUC = kL
        
    """
    Function: normalize the thermal conductivity values
    """
        
    def normalize(self):
         kL = self.kL 
         return kL/(max(kL))
    
    
    """
    Function: Born von Karmann dispersion
    """
    def bvk_kappa(self):
        vs = self.vs
        vol = self.vol
        comp = self.data.comp
        c_len = 100
        kL = np.zeros(c_len)
        kmax = (6*pi**2/vol)**(1/3) #same integration limit.. in terms of k
        #print(kmax)
        def vg(k): 
            vg = vs*cos(pi/2* k/kmax)
            return vg
        def vp(k): 
            vp = vs*(2/pi)*(kmax/k)*sin(pi/2*k/kmax)
            return vp
        def vg3(k):
            vg3 = (vs*(cos((pi/2)* (k/kmax))))**3
            return vg3
       # print(vg3(kmax))
        def omega(k): 
            omega = vp(k)*k
            return omega
        vg3_avg, other = integrate.quad(vg3, 0, kmax)
#        omegaD = (6*pi**2/vol)**(1/3)*vs
#        omegaRange = np.linspace(0,omegaD, 100)
#        plt.plot(omegaRange, )
        #print(vg3_avg)
        #kappa0 = kb/(2*pi**2)*vg2_avg # I have kappa0 though...
        crange = np.linspace(comp[0], comp[len(comp) -1], c_len)
        i = 0
        for c in crange:
            self.bfz(c)
            kap_pure = self.data.kappa[0]*(1-c) + self.data.kappa[len(self.data.kappa)-1]*c
            a = vol/(4*pi)*self.gamma
            b = kb/(2*pi**2*kap_pure)*vg3_avg
            if i==0:
                kL[i] = self.data.kappa[0]
            else:
                def kInt(k):
                    kInt = vg(k)**3*(1/(1+ a*(omega(k)**2/b)))
                    return kInt
                kappaInt, other = integrate.quad(kInt, 0, kmax)
                print(kappaInt)
                kL[i] = kap_pure*(1/vg3_avg)*kappaInt
            i=i+1
        self.kL = kL