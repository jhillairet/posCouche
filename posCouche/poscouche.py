# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
import warnings

from pywed import tsbase, tsmat, PyWEDException

# Hack to be able to run this fil directly
if __name__ == '__main__':
    from formeTS import vacuum_vessel, LCFS
    from ic_utils import IC_resonance_radius_ripple, IC_resonance_radius, IC_resonance_frequency, WEST_toroidal_field
else:
    from . formeTS import vacuum_vessel, LCFS
    from . ic_utils import IC_resonance_radius_ripple, IC_resonance_radius, IC_resonance_frequency, WEST_toroidal_field

    
from scipy.constants import m_p, e

def poscouche(shot=None, time=None, Itor=None, freq=55.0, n=1, species='H'):
    """   
    Display the locations of Ion Cyclotron Resonance 
    layers for a given minority ion species for Tore Supra/WEST.
    
    Arguments:
    - shot: shot number (default:None)
    - time: shot time [s] (default: time at the middle of the shot)  
    - Itor: current in toroidal coils [A] (default: 1250)
    - freq: RF frequency [MHz] (default: 55)
    - n : harmonic number (1,2,3,...) (default: 1)
    - species : 'H', 'D', 'T', '3He', '4He' (default: 'H')
    
    Returns nothing
    
    Examples:
       >> poscouche(freq=60)
       >> poscouche(Itor=900)
       >> poscouche(Itor=900, freq=60)
       >> poscouche(Itor=900, freq=60, n=2)
       or:
       >> poscouche(shot=47979)
       >> poscouche(shot=47979, time=10)
    
    @author: Vincent Basiuk, Julien Hillairet
    """      
    # NB: WEST configuration > 48451
    
    # large and minor radius in meter
    R0 = 2.37
    a = 0.80
    
    # set default values if shot number not passed
    if shot is None:
        # samin: small radius
        xa = 0.8
        ta = 1
        # srmaj: large radius
        x = 2.4
        t = 1
        #axe = R0 + d
        # sd0mag
        d = 0.1
        t1 = 1
        # Antenna frequencies
        if not isinstance(freq, list):
            freq = [freq, freq, freq]
        ant = [0, 0, 1]
        # Antenna positions
        pos = [3.3, 3.3, 3.3]     
        # shot time 
        time = 5    
        
        if Itor is None:
            Itor = 1250 # A

    else:        
        try:
            xa,ta = tsbase(shot, 'samin',  nargout=2)
            x,t = tsbase(shot, 'srmaj', nargout=2)
            d,t1 = tsbase(shot, 'sd0mag', nargout=2)
            freq, ant = tsbase(shot, 'SFREQFCI', nargout=2)
            pfci, t_pfci = tsbase(shot, 'SPUISS', nargout=2)
            pos = tsbase(shot, 'SPOSFCI', nargout=1)
            
            # Set the shot time to the time of max ICRH power
            # Otherwise set the it to the middle of the shot
            if time is None:
                if pfci is not None:
                    time = t_pfci[np.where(pfci == np.max(pfci))]
                else: 
                    if shot < 28540:
                        time = np.max(t)/2
                    else:
                        time = np.max(ta[xa>0])/2

            # Plasma properties at time
            R0 = x[np.argmin(np.abs(t - time))]
            a = xa[np.argmin(np.abs(t - time))]
            #axe = R0 + d[np.argmin(np.abs(t - time))]

            # Current in toroidal coils at the time
            if shot < 20100:
                Itor = tsmat(shot, 'EXP=T=S;EXP=T=S;ITOR')
            if shot < 20600:
                Itor = tsmat(shot, 'EXP=T=S;GENERAL;ITOR')
            else:
                Itor, tsi = tsbase(shot, 'SITOR', nargout=2)
            Itor = Itor[np.argmin(np.abs(tsi - time))] 
                        
        except PyWEDException:
            raise PyWEDException('Error with tsbase')
    
    fig, ax = plt.subplots()
    
    # Plot vacuum chamber and plasma LCFS profiles.
    # wall profile available only for shot > 28540        
    if shot is not None:
        R_vv, Z_vv = vacuum_vessel(shot)
        ax.plot(R_vv, Z_vv, 'k', lw=2)

        try: 
            Rext, Zext, t = LCFS(shot)
            ind_t = np.argmin(np.abs(t - time))
            ax.plot(Rext[ind_t,:], Zext[ind_t,:], 'b', lw=1)
        except ValueError:
            warnings.warn('No LCFS data available, passing.')
    else:
        print('No shot number provided: using WEST vacuum vessel')
        R_vv, Z_vv = vacuum_vessel(shot=50000)
        ax.plot(R_vv, Z_vv, 'k', lw=2)       
        
    ax.axis('equal')
    ax.set_xlabel('R [m]', fontsize=14)
    ax.set_ylabel('Z [m]', fontsize=14)
    ax.set_title('Shot #{}@t={}s'.format(shot, time), fontsize=14)
    
#    # plot antenna positions and frequencies
#    ax.plot([pos[0], pos[0]],[-0.6, 0.6], 'c',
#             [pos[0], pos[0]+0.1], [-0.6, -0.6],'c',
#             [pos[0], pos[0]+0.1], [0.6, 0.6],'c', lw=2)
#    ax.plot([pos[1], pos[1]], [-0.6, 0.6],'y',
#             [pos[1], pos[1]+0.1], [-0.6, -0.6],'y',
#             [pos[1], pos[1]+0.1], [0.6, 0.6],'y', lw=2)
#    ax.plot([pos[2], pos[2]], [-0.6, 0.6],'c--', 
#             [pos[2], pos[2]+0.1], [-0.6, -0.6],'c--', 
#             [pos[2], pos[2]+0.1], [0.6, 0.6],'c--', lw=2)
#             
#    ax.text(x=pos[0]+0.1, y=0.5, s='Q1: {} m\n   {} MHz'.format(
#            int(pos[0]*100)/100, int(freq[0]*100)/100))
#    ax.text(x=pos[1]+0.1, y=0.3, s='Q4: {} m\n   {} MHz'.format(
#            int(pos[1]*100)/100, int(freq[1]*100)/100))
#    ax.text(x=pos[2]+0.1, y=0.1, s='Q5: {} m\n   {} MHz'.format(
#            int(pos[2]*100)/100, int(freq[2]*100)/100))  
    
    # domain space 
    R = np.linspace(1.5, 3.5, 501) # arbitrary values for min and max. 
    z = np.linspace(-a, a, 101)
    RR, ZZ = np.meshgrid(R, z)   
    
    B = WEST_toroidal_field(Itor, R)

#    fp = e*B/(2*np.pi*m_p*1e6)
    fp = IC_resonance_frequency(B, species)    
    # Plot the resonance layers for the three harmonics k=1,2,3
    colours = iter(cm.copper(np.linspace(0,1,len(freq))))
    
    # only perform the loop for frequencies which are different 
    freq = np.unique(freq) 
    for f in freq:
        # index at which the resonance condition 
        # starts being satisfied
        ind = np.argwhere(f/n >= fp)
        if ind.size:
            res_layerH = R[ind[0]]

            R_wo_ripple, R_ripple_max, R_ripple_min = IC_resonance_radius_ripple(
                              RR, ZZ, Itor, f, n, species=species)    

            cur_col = next(colours)
                
            ax.plot(R_ripple_max, z, color=cur_col)   
            ax.plot(R_ripple_min, z, color=cur_col)   
                
            print('Resonance radius for {} = {} m @{} MHz (n={})'.format(species, R_wo_ripple, f, n))
            ax.axvline(R_wo_ripple, color=cur_col, ls='--', lw=2)
            ax.text(x=R_wo_ripple+0.05  , y=0, s='{}@{} MHz (n={})'.format(species,f,n), 
                    rotation=90, fontsize=14, color='b')
            ax.fill_betweenx(z, R_ripple_min, R_ripple_max, alpha=0.2, color=cur_col)
    fig.show()

if __name__ == '__main__':
    #Possible usage of the posCouche function :
    #poscouche(Itor=1250, freq=60, n=1, species='H')
    #poscouche(freq=55, Itor=1250)
    #poscouche(freq=[50, 55, 60], Itor=1250)
    #poscouche(shot=47979, time=10)       
    Itor = input('Entrez la valeur du courant toroidal [A] (default: 1250 A) :')
    if not Itor:
        Itor = 1250
    else: 
        Itor = float(Itor)
    
    freq = input('Entrez la valeur de la frequence HF [MHz] (default: 55.5 MHz) :')
    if not freq:
        freq = 55.5
    else:
        freq = float(freq)

    print('Espèce minoritaire : H. Harmonique par defaut: n=1.')
    
    poscouche(Itor=float(Itor), freq=float(freq), n=1, species='H')
    
    

    
    
