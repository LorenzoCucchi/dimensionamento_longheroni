# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 14:11:22 2021

@author: Lorenzo Cucchi


"""

import matplotlib.pyplot as plt
import numpy as np


def area_calculator(wing):
    areasec = []
    span = 0
    for x in range(1,len(wing),1):
       a = wing[x-1]
       b = wing[x]
       areasec.append((a[0]+b[0])*b[1]/2);
       span = span+b[1]
    area = sum(areasec)*2/(10**6)
    span = (span+wing[0][1])*2    
    return area,span
       
def spar_width(wing):
    ds = 1
    b = []
    s = []
    for x in range(1,len(wing)):
        b1 = wing[x-1][3]
        b2 = wing[x][3]
        L = wing[x][1]
        s1 = wing[x-1][2]*wing[x-1][0]/100
        s2 = wing[x][2]*wing[x][0]/100
        tangenteb = (b1-b2)/(2*L)
        tangentes = (s1-s2)/(2*L)
        for i in range(0,wing[x][1],ds):
           b.append((L-i)*tangenteb*2 + b2)
           s.append((L-i)*tangentes*2 + s2)
           
    return b,s  


def sforzi(q,span,wing,sigma_comp,sigma_tens,tau,nn):
    Mb = []
    Q = []
    Wb_comp = []
    Wb_tens = []
    F = []
    L = int(span/2)
    Li = wing[0][1]
    for x in range(Li,L):
        Mb.append((q*(L - x)**2)/2)
        Wb_comp.append((q*(L-x)**2)*(1-nn)/sigma_comp)
        Wb_tens.append((q*(L-x)**2)*nn/sigma_tens)
        Q.append(q*(L-x))
        F.append(q*(L-x)/tau)
    return Mb,Q,Wb_comp,Wb_tens,F 


def dimensions(b,s,Wb_comp,F,rovArea,nn):
    
    t1 = []
    Rov_upp = []
    Rov_low = []
    H_upp = []
    H_low = []
    for x in range(len(b)):
        hc=np.arange(0,0.2,0.0001)
        ht=nn-np.sqrt(nn**2-2*hc+hc**2+2*nn*hc)
        W_prov = b[x]*((hc**3+ht**3)/12 + hc*(s[x]*(1-nn - hc/2))**2+ ht*(s[x]*(nn-ht/2))**2)/(1-nn)
        ind=np.abs(W_prov - Wb_comp[x]).argmin()
        
        H_upp.append(s[x]*np.take(hc,ind))
        H_low.append(s[x]*np.take(ht,ind))  
    
        Rov_upp.append(H_upp[x]*b[x]/rovArea)
        Rov_low.append(H_low[x]*b[x]/rovArea)
        
        t1.append(F[x]/(2*(b[x]-H_upp[x]-H_low[x])))
    return H_upp,H_low,t1,Rov_low,Rov_upp,ind  

        
        
        
sections = []
# units are in millimeters

sections.append([255.,  81, 12., 50.]) #50
sections.append([240., 495, 11, 50.]) #50
sections.append([211., 905, 8.7, 40.]) #40
sections.append([163., 705, 8.7, 30.]) #30
sections.append([107., 560, 8.7, 20.])#20
sections.append([ 76., 196, 8.7, 15.])#15





area,span = area_calculator(sections)

" massa in kg escluse le ali "
mass = 9.


# load factor in g
n = 20
g = 9.8
weight = mass*g*n

" dati materiale N/mm2 "
sigma_tens = 750
sigma_comp = 500
sig_frac = sigma_comp/sigma_tens
tau = 10
rovArea = 2
nn = 1/(sig_frac+1)

" reazioni vincolari "
Mbmax = weight*span/12
q = 12*Mbmax/span**2
Qmax = q*span/2

" caratterizzazione larghezza longherone "
b,s = spar_width(sections)
  
" calcolo sforzi "
Mb,Q,Wb_comp,Wb_tens,F = sforzi(q, span, sections,sigma_comp,sigma_tens,tau,nn)

" calcolo dimensioni longherone e rovings"
H_upp,H_low,t1,Rov_low,Rov_upp,ind = dimensions(b,s,Wb_comp,F,rovArea,nn)



ROV_upp=[]
ROV_low=[]
for x in range(len(b)):
   ROV_upp.append(np.ceil(Rov_upp[x]))
   ROV_low.append(np.ceil(Rov_low[x]))



fig1, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
ax1.plot(Rov_low)
ax1.plot(Rov_upp)
ax1.plot(ROV_low)
ax1.plot(ROV_upp)
ax1.set_title("Roving distribution and thickness")
ax1.set_ylabel("Number of roving")
ax1.legend(['Lower Roving','Upper roving'])

ax2.plot(H_low)
ax2.plot(H_upp)
ax2.set_xlabel("Span position")
ax2.set_ylabel("Thickness")
ax2.legend(['Lower Thicknes','Upper Thickness'])

plt.savefig('rovings.svg')