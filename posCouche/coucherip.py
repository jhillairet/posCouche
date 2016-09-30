# -*- coding: utf-8 -*-
# NumPy
import numpy as np
    
def coucherip(R, Z, B0, freq, n, ep, A):
    """
    Calculates the radius of the Ion Cyclotron Resonance layer, 
    taking into account the magnetic ripple. 
    Three values are returned: 
     - the resonance condition  
     - the maximum radius (under in-between coils)    
     - the minimum radius (under the coils) 
    
    The ripple function is calculated from an analytical 
    from V.Basiuk et al., Fusion Technology 26 (Nov 1994) p.222-226
    
    Arguments:
    - R: Large radius [m]
    - Z: vertical position [m]
    - B0: magnetic field at plasma center [T]
    - freq: RF frequency [MHz]
    - n: harmonic number (1, 2, ...)
    - ep: +1 or -1
    - A: mass number, ie. total number of protons and neutrons
    
    Returns:
    - res_cond: the resonance condition f_ci - f_rf/n with ripple [Hz]
    - R_ripple: Radius of the resonance layer [m]
    - R_wo_ripple: Radius of the resonance layer without ripple [m]

    
    Authors: V.Basiuk, J.Hillairet
    """
    # Import some physical constants from scipy.constants
    from scipy.constants import m_p, e, pi 
    # Convert into numpy array, because we need some array methofs
    R = np.array(R)
    Z = np.array(Z)  
    
    R0 = 2.37 # Tore Supra Major Radius
    B = B0*R0/R # Magnetic field at R
    fci = e*B/(2*pi*m_p*A) # cyclotron frequency at R
    
    X = R - 2.04
    Y = 0.52*X + 1
    b = 0.26
    b2 = 2*b**2
    b4 = 4*b**2
    rc = np.sqrt(np.abs(Y - np.sqrt(np.abs(Y**2 - b4*(X**2+Z**2))))/b2)
    rip = 1 + ep*2.2e-4*np.exp(rc*(5 + rc*1.6))
      
    # The resonance layers correspond to radius 
    # which minimize the conditions:
    res_cond = fci*rip - freq*1e6/n 
    res_cond_wo_rip = fci - freq*1e6/n      

    # TODO ? check if a solution exists
          
    # Radii which minimize the above conditions
    # flat is necessary is the dimension of the R array > 1
    # axis=1 leads to return the radii for each Z values
    if np.ndim(res_cond) == 1:
        R_ripple = R.flat[np.argmin(np.abs(res_cond))]
        R_wo_ripple = R.flat[np.argmin(np.abs(res_cond_wo_rip))]
    else:
        R_ripple = R.flat[np.argmin(np.abs(res_cond), axis=1)]
        R_wo_ripple = R.flat[np.argmin(np.abs(res_cond_wo_rip), axis=1)]
    
    return res_cond, R_ripple, R_wo_ripple
    

# The following code is run if one executes this file directly      
if __name__ == '__main__':
    from matplotlib.pyplot import *
    # Generate a R,Z grid
    z = np.linspace(-1.2, 0, 50)    
    R, ZZ = np.meshgrid(np.linspace(2.0, 2.9, 501), z)
    res_cond, R_ripple, R_wo_ripple= coucherip(R, ZZ, 3.785, 57, 1, -1, 1)
    print(R_ripple, R_wo_ripple)
    figure()
    plot(R_ripple, z)
    plot(R_wo_ripple, z)
    