# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import argrelmin
from scipy.constants import m_p, m_n, e, pi

def WEST_toroidal_field(Itor=1250, R=2.37):
    """
    Returns the WEST toroidal field magnitude as a function of the radius R
    and the current flowing in the toroidal coils
    
    Arguments:
        - Itor: current in the toroidal coils [A] (default: 1250)
        - R: radius [m] (default: 2.37)
        
    Returns:
        - B: Magnetic field at R [T]
    """
    return 0.0073*Itor/R

def ion_mass_and_charge(species='H'):
    """
    Returns the fully ion isotope species mass in [kg]
    
    m, q = ion_mass_and_charge(species)
    
    Argument:
        - species : 'H', (or '1H', 'hydrogen', 'proton', 'protium'), 
                    'D', (or '2D', 'deuterium', 'deuteron')
                    'T', (or '3H', 'tritium')
                    'He', (or '4He', 'helium')
                    '3He', (or 'helion')
    Returns:
        - m: ion mass [kg]
        - q: ion electric charge [C]
    """
    SPECIES=str.upper(species)
    if SPECIES in ('H', '1H', 'HYDROGEN', 'PROTIUM'):
        A = 1; Z = 1
    elif SPECIES in ('D', '2H', 'DEUTERIUM', 'DEUTERON'):
        A = 2; Z = 1
    elif SPECIES in ('T', '3H', 'TRITIUM'):
        A = 3; Z = 1
    elif SPECIES in ('HE', '4HE', 'HELIUM'):
        A = 4; Z = 2
    elif SPECIES in ('3HE', 'HELION'):
        A = 3; Z = 2
    else:
        raise ValueError('Incorrect species argument: {}'.format(species))
        
    m = Z*m_p + (A-Z)*m_n
    q = Z*e
    return m, q

def IC_resonance_frequency(B=3.86, species='H'):
    """
    Returns the fundamental cyclotron resonance frequency (in MHz)
    
    Arguments:
        B: magnetic field magnitude [T] (default: 3.86)
        species: '1H', '2H', '3H', '4He', '3He' (more possibilities, cf ion_mass_and_charge definition)
        
    Returns:
        f: fundamental cyclotron resonance frequency [MHz]
    """
    m, q = ion_mass_and_charge(species)
    return q*B/(2*pi*m)/1e6

def IC_resonance_radius(Itor=1250, f=55, n=1, species='H'):
    """
    Calculates the radius of the Ion Cyclotron resonance layer.
    
    Arguments:
        Itor: current in the toroidal coils [A] (default: 1250)
        f: RF frequency [MHz] (default: 55)
        n: harmonic number (1, 2, ...) (default: 1)
        species: '1H', '2H', '3H', '4He', '3He' (more possibilities, cf ion_mass_and_charge definition)
    
    Returns:
        R_ic: Ion Cyclotron Resonance radius [m]
        
    """
    # Toroidal field at R0
    R0 = 2.37
    B0 = WEST_toroidal_field(Itor=Itor, R=R0)
    # ion mass 
    m, q = ion_mass_and_charge(species)
    # cyclotron resonance radius
    R_ic = n*(q/m)*(B0*R0)/(2*pi*f*1e6) 
    
    return R_ic
    
def IC_resonance_radius_ripple(R, Z, Itor=1250, freq=55, n=1, species='H'):
    """
    Calculates the radius of the Ion Cyclotron Resonance layer, 
    taking into account the magnetic ripple.
    
    Arguments:
    - R: Large radius [m]
    - Z: vertical position [m]
    - Itor: current in the toroidal coils [A] (default: 1250)
    - freq: RF frequency [MHz]
    - n: harmonic number (1, 2, ...)
    - species: '1H', '2H', '3H', '4He', '3He' (more possibilities, cf ion_mass_and_charge definition)
    
    Returns:
    - R_wo_ripple : radius of the resonancer layer without ripple [m]    
    - R_ripple_min: minimum radius of the resonance layer (under the coil) [m]
    - R_ripple_max: maximum radius of the resonance layer (under in-between coils) [m]

    References: The ripple function is calculated from an analytical 
    from V.Basiuk et al., Fusion Technology 26 (Nov 1994) p.222-226
     
    Authors: V.Basiuk, J.Hillairet
    """
    # Convert into numpy array, because we need some array methods
    R = np.asarray(R)
    Z = np.asarray(Z)  
    

    B = WEST_toroidal_field(Itor, R)
    B0 = WEST_toroidal_field(Itor)
    
    m, q = ion_mass_and_charge(species)
    
    R_wo_ripple = IC_resonance_radius(Itor, freq, n, species)
    
    fci = q*B/(2*pi*m) # cyclotron frequency at R
    
    X = R - 2.04
    Y = 0.52*X + 1
    b = 0.26
    b2 = 2*b**2
    b4 = 4*b**2
    rc = np.sqrt(np.abs(Y - np.sqrt(np.abs(Y**2 - b4*(X**2+Z**2))))/b2)
    rip_min = 1 + (-1)*2.2e-4*np.exp(rc*(5 + rc*1.6))
    rip_max = 1 + (+1)*2.2e-4*np.exp(rc*(5 + rc*1.6))
    
    # The resonance layers correspond to radius 
    # which minimize the conditions:
    res_conds = {'min': fci*rip_min - freq*1e6/n, 
                 'max': fci*rip_max - freq*1e6/n} 
    res_cond_wo_rip = fci - freq*1e6/n      

    # output array initialisation
    R_ripple = np.zeros(Z.shape[0])
    R_wo_ripple_ = np.zeros(Z.shape[0])
    
    # find if a solution exists Z by Z (ie. line by line)
    # One could have made the calculation directly for the full array,
    # but then you may find incorrect solutions (for high Z) 
    # because sometime there is simply no solution at all
    for cond in res_conds:
        for idz,z in enumerate(Z):   
            res_cond_z = res_conds[cond][idz,:]
    
            # Radii which minimize the resonance condition wo ripple
            # (Gets the index of the closest value to 0)
            R_wo_ripple_[idz] = R_wo_ripple# R.flat[np.argmin(np.abs(res_cond_wo_rip))]
            
            # Check if a solution exists for the ripple case.
            # If a continuous function has values of opposite sign inside an interval,
            # then it has a root in that interval (Bolzano's theorem)                        
            if (np.sign(np.min(res_cond_z)) == -1) & (np.sign(np.max(res_cond_z)) == +1):
                # Depending of the radius R range, you may find few different solutions 
                # (which are in fact harmonics resonance layers)
                # The solution we look for is the one closest to the resonance layer wo ripple
               
                # Get the relative minimas of the resonance condition
                idx_relmin, = argrelmin(np.abs(res_cond_z)) # returns a tuple of ndarray
                # Select the index which corresponding values is the closest 
                # of resonance layer wo ripple
                try:
                    idx_res = idx_relmin[np.argmin(np.abs(R.flat[idx_relmin] - R_wo_ripple_[idz]))]
                    R_ripple[idz] = R.flat[idx_res]
                except ValueError:
                    # no resonance condition satisfied in the domain range
                    R_ripple[idz] = np.NAN
            else:
                # no resonance condition satisfied in the domain range
                R_ripple[idz] = np.NAN

        if cond is 'min':
            R_ripple_min = R_ripple.copy()
        elif cond is 'max':
            R_ripple_max = R_ripple.copy()

    return R_wo_ripple, R_ripple_min, R_ripple_max
    
# The following code is run if one executes this file directly      
if __name__ == '__main__':
    from matplotlib.pyplot import *
    # Generate a R,Z grid
    z = np.linspace(-1.2, 1.2, 101)
    r = np.linspace(1.5, 3.5, 501)     
    R, ZZ = np.meshgrid(r, z)
    
    Itor = 1250
    freq=55
    ns=[1,2,3]
    species = 'D'
    
    figure(1)
    clf()
    
    for n in ns:    
        R_wo_ripple, R_ripple_min, R_ripple_max = IC_resonance_radius_ripple(
                              R, ZZ, Itor, freq, n, species)
        
        plot(R_ripple_max, z)
        plot(R_ripple_min, z)
        axvline(R_wo_ripple, ls='--', lw=2)
     
    axis('equal')
