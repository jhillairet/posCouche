# -*- coding: utf-8 -*-
"""
Plot the Ion Cyclotron resonance layer versus Frequency 
"""
from scipy.constants import pi, m_p, e
import numpy as np
from matplotlib.pyplot import *
matplotlib.rcParams.update({'font.size': 16})
from coucherip import coucherip

def resonance_R (f, R0, B0, A, n):
    '''
    Return the cyclotron resonance frequencies in [Hz] 
    at LFS, R0 and HFS for a given:
        - f: RF frequency in [Hz] 
        - R0: large radius (axis) [m]
        - B0: toroidal magnetic field at R0
        - A: Atomic mass number
        - n: resonance index        
    '''
    R_c = (n*e*B0*R0)/(f*2*pi*m_p*A)
    return R_c

# Atomic mass
A_H = 1
A_D = 2
A_He3= 3

## Tore Supra parameters    
#R0 = 2.4 # m
#a = 0.72 # m
#source_frequencies = [42, 48, 57, 63, 76]

# WEST parameters    
#R0 = 2.5 # m
a = 0.5 # m
B0 = 3.7 # T 
source_frequencies = [48, 53, 55.5, 57, 63]

# Calculates the cyclotron frequencies vs f
f_RF = np.linspace(40e6, 65e6, 30)

R_cis_2H, R_cis_1H, R_cis_He3, R_cis_2He3 = [], [], [], []
for f in f_RF:
    R_cis_1H.append(resonance_R(f, R0, B0, A=A_H, n=1))
    R_cis_2H.append(resonance_R(f, R0, B0, A=A_H, n=2))
    R_cis_He3.append(resonance_R(f, R0, B0, A=A_He3, n=2))  
    R_cis_2He3.append(resonance_R(f, R0, B0, A=A_He3, n=4))  
R_cis_1H = np.array(R_cis_1H)
R_cis_2H = np.array(R_cis_2H)
R_cis_He3 = np.array(R_cis_He3)
R_cis_2He3 = np.array(R_cis_2He3)

# Plotting
fig, ax=subplots(1,1)
# 1H
#fill_between(f_RF, R_cis_1H[:,0], R_cis_1H[:,2], alpha=0.2, color='b')
plot(f_RF/1e6, R_cis_1H, lw=2, color='b')
text(45,2.8, '$H$', color='b', fontsize=30)

# 2H
#fill_between(f_RF, R_cis_2H[:,0], R_cis_2H[:,2], alpha=0.2, color='g')
plot(f_RF/1e6, R_cis_2H, lw=2, color='g')
text(2,70, '$2H$', color='g', fontsize=30)

# 1He3
#fill_between(f_RF, R_cis_2H[:,0], R_cis_2H[:,2], alpha=0.2, color='g')
plot(f_RF/1e6, R_cis_He3, lw=2, color='r')
text(45,2.2, '$^3He$', color='r', fontsize=30)

# 2He3
#fill_between(f_RF, R_cis_2He3[:,0], R_cis_2He3[:,2], alpha=0.2, color='r')
plot(f_RF/1e6, R_cis_2He3, lw=2, color='c')
text(60,3.2, '$2^3He$', color='c', fontsize=30)

ax.set_xlabel('$f_{RF}$ [MHz]', fontsize=18)
#ax.set_xticks([1, 2, 3, 4, 5])
ax.set_ylabel('R [m]', fontsize=18)
ax.grid(True)

ax.set_ylim(1.6, 3.4)
plt.axhspan(1.6, 2, color='gray', alpha=0.25)
ax.hlines(2.0, 40, 65, color='k', linestyles='--', lw=2)
ax.hlines(2.5, 40, 65, color='gray', linestyles='--', lw=2)
ax.hlines(3.0, 40, 65, color='k', linestyles='--', lw=2)
plt.axhspan(3, 3.4, color='gray', alpha=0.25)

title('$R_0$=2.5 m, $B_0$=3.7 T')
#R_1H = []
R_ripple = []
for f in f_RF:
    res_cond, R_ripple1, R_wo_ripple = coucherip(np.array([R0]), np.array([0]), B0, f/1e6, 1, +1, A_H)
    R_ripple.append(R_ripple1)

for f_s in source_frequencies:
    ax.vlines(f_s, 1.6, 3.4, linestyle=':', color='navy')
#    plt.axvspan(f_s-2, f_s+2, color='b', alpha=0.25)
#    res_cond, R_ripple2, R_wo_ripple = coucherip(np.array([R0]), np.array([0.5]), b, 55, 1, -1, A_H)
#    R_1H.append([R_ripple1, R_ripple2])
#    
#R_1H = np.array(R_1H)