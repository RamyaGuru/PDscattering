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
        self.comp = data[:,0]
        self.kappa = (data[:,1])
    
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
            #self.data = kLcurve(fileID)
    
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
            gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
            gamma = gamma + gamma_s
            print("For site %d, the gamma is: %4.8f" % (n, gamma_s))
        self.gamma = gamma
            
            
#    def bfz(self, c):
#        stoich = self.stoich
#        subst = self.subst
#        mass = self.mass
#        #m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
#        natoms = sum(stoich)
#        gamma = 0
#        for n in range(len(mass)):
#            #msite = subst[n]*c + mass[n]*(1-c)
#            gamma_s = 0
#            for i in range(len(subst[n])): 
#                msite = subst[n][i]*c + mass[n]*(1-c)
#                m_uc = np.dot(stoich, sum(subst[n])/len(subst[n]))*c + np.dot(stoich, mass)*(1-c)
#                gamma_s = gamma_s + (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n][i])/msite)**2
#                print(gamma_s)
#            gamma = gamma + gamma_s
#            self.gamma = gamma_s

    
    """
    Function: mass fluctuation parameter calculation with the PUC parameter 
    """
 #Need to adjust: still treating sublattice independently    
    def puc(self, c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass       
        gammaPUC = 0
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c)
        #print(m_uc)
        natoms = sum(stoich)
        p_full = 0
        for n in range(len(mass)):
            gamma_s = 0
            for m in range(0,stoich[n]+1):
                #msite = subst[n]*c + mass[n]
                m_uc_nm = subst[n]*m + mass[n]*(stoich[n] - m) + np.dot(np.delete(stoich, n), np.delete(mass, n))*(1-c) + np.dot(np.delete(stoich, n), np.delete(subst, n))*(c)
               # print(m_uc_nm)
                p = binom.pmf(m, stoich[n], c)
                gamma_nm = binom.pmf(m, stoich[n], c)*(1-(m_uc_nm/m_uc))**2
                #print(binom.pmf(m, stoich[n], c))
                gamma_s = gamma_s + gamma_nm
                print("Concetration: %2.2f, site number: %2d, site Gamma: %4.8f" % (c, m, gamma_s))
                p_full = p_full + p
            gammaPUC = (gammaPUC+ gamma_s) 
        print('full p: %4.8f' % p_full)
        self.gamma = gammaPUC*natoms
    
#    def puc_multi(self, c):
#        stoich = self.stoich
#        subst = self.subst
#        mass = self.mass       
#        gammaPUC = 0
#        m_uc = np.dot(stoich, sum(subst[n])/len(subst[n]))*c + np.dot(stoich, mass)*(1-c)
#        natoms = sum(stoich)
#        for m in range(len(mass)):
#            for n in range(len)

#So far, this will only work for an A and B lattice-- try to extend to more complex crystal 
# Re-write as a recursive function


    def puc_proper2(self,c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass
        gammaPUC = 0
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
        count = [0]
        count2 = [0]
        func = [0]*sblatt
        loop = [0]*sblatt
        loop2 = [1]*sblatt
        #Have a counter for which function call I am in
        def recloop(a):
            func[a] +=1
            print(func)
            for m in range(stoich[a]+1):
                loop[a]+=1
                loop2[a] =0
                print(loop)
                print('loop2')
                print(loop2)
                print('m_def  and m here: %4.8f %d' % (m_def[0], m))
                pUC_old = pUC[0] #this shouldn't get called for every m
                m_def_old = m_def[0]
                pUC[0] = pUC[0]*binom.pmf(m,stoich[a],c)
                m_def[0] = m_def[0] + subst[a]*m +mass[a]*(stoich[a] - m)
                if a+1<sblatt:
                    #print(func)
                    #print('a is %d' % (a+1))
                    pUC_pregen[a+1] = pUC[0]
                    m_def_pregen[a+1] = m_def[0]
                    print(m_def_pregen)
                    loop2[:] = [1]*sblatt
                    recloop(a+1)
                else:
                    #print(count[0])
                    #count[0] +=1
                    gamma_nm = (1-(m_def[0])/m_uc)**2
                    print(m_def)
                    print(gamma_nm)
                    print(pUC[0])
                    GammaTot[0] = GammaTot[0] + pUC[0]*gamma_nm
                    #Have this if statement before the calculation?
                    if m==stoich[a]: #NEED TO MODIFY: so that this will work for more than 2 sublattices... needs to know whether the loop is on the last sublattice
                        print('a = %d' % a)
                        pUC[0] = pUC_pregen[a-1] #should subtract based on the function call number... 
                        m_def[0] = m_def_pregen[a-1] 
                    else:
                        pUC[0] = pUC_old
                        m_def[0] = m_def_old
                        print(count2[0])
                        count2[0] +=1
                    loop2[:] = [1]*sblatt
        recloop(0)
        print(GammaTot[0])
        self.gamma = GammaTot[0]*natoms
                        
          
    def puc_proper(self,c):
        stoich = self.stoich
        subst = self.subst
        mass = self.mass       
        gammaPUC = 0
        m_uc = np.dot(stoich, subst)*c + np.dot(stoich, mass)*(1-c) 
        natoms = sum(stoich)
        full_p = 0
        for m in range(stoich[0] +1):
            pA = binom.pmf(m, stoich[0], c)
            for n in range(stoich[1] + 1):
                pB = binom.pmf(n, stoich[1], c)
                ptot = pA* pB
                m_uc_nm = subst[0]*m + mass[0]*(stoich[0] -m) + subst[1]*n + mass[1]*(stoich[1]-n)
                gamma_nm = ptot*(1-(m_uc_nm/m_uc))**2
                gammaPUC = gammaPUC + gamma_nm
                print(ptot)
                print('The gamma = %4.8f for n = %d and m = %d and m_uc_nm = %4.8f and m_uc = %4.8f' % (gamma_nm, n, m, m_uc_nm, m_uc))
                full_p = full_p + ptot
        print('full_p: %4.8f' % full_p)
        self.gamma = gammaPUC*natoms
        
        
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
                gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*(-(abs(subst[n] - mass[n]))/msite - 2)**2
            else:
                msite = subst[n]*c + mass[n]*(1-c)
                gamma_s = (stoich[n]/natoms)*c*(1-c)*(msite/(m_uc/natoms))**2*((mass[n] - subst[n])/msite)**2
            gamma = gamma + gamma_s
        self.gamma = gamma


    """
    Function: calculate the thermal conductivity from gamma
    """
    def KappaL(self): #for u-curves... considers both end-members
        num = sum(self.stoich)
        comp = self.data.comp
        vol = self.vol
        vs = self.vs
        prefix = (6**(1/3)/2)*(pi**(5/3)/kb)*(vol**(2/3)/vs)
        c_len = 100
        kL = np.zeros(c_len)
        crange = np.linspace(comp[0], comp[len(comp) - 1], 100)
        i = 0
        for c in crange:
            self.bfz_vac(c)
            kap_pure = self.data.kappa[0]*(1-c) + self.data.kappa[len(self.data.kappa)-1]*c
            u = (prefix*self.gamma*kap_pure)**(1/2)
            if i==0:
                kL[i] = self.data.kappa[0]
            else:
                kL[i] = kap_pure*np.arctan(u)/u
            i = i+1
        self.kL = kL
        
    """
    Function: calculate the thermal conductivity when both the 
    """
        
    def normalize(self):
         kL = self.kL 
         return kL/(max(kL))
    