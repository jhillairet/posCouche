# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import warnings

from pywed import tsbase, tsmat, PyWEDException

from formeTS import formeTS
from ic_utils import IC_resonance_radius_ripple

from scipy.constants import m_p, e

def poscouche(shot=None, time=None, Itor=None, freq=55):
    """
    Poscouche
    
    Calculates and display the locations of Ion Cyclotron Resonance 
    layers for Hydrogen and Deuterium.
    
    Keywords Arguments:
    - shot: shot number
    - time: shot time [s] (default: time at the middle of the shot)  
    - Itor: current in toroidal coils [A] (default: 1250)
    - freq: RF frequency [MHz] (default: 55)
    
    Returns nothing
    
    Examples:
       >> poscouche(freq=60)
       >> poscouche(Itor=900)
       >> poscouche(Itor=900, freq=60)
       >>
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
    
    fig, ax = plt.subplots(1,1)
    
    # Plot vacuum chamber and plasma LCFS profiles.
    # wall profile available only for shot > 28540        
    if shot is not None:
        if shot > 28540:
            Rparoi, Zparoi, Rext, Zext, t = formeTS(shot)
            ind_t = np.argmin(np.abs(t - time))
            ax.plot(Rparoi, Zparoi, 'k', lw=2)
            ax.plot(Rext[ind_t,:], Zext[ind_t,:], 'b', lw=1)
    else:
        warnings.warn('No vacuum chamber information: ideal wall profile plotted')
        theta = np.linspace(0, 2*np.pi, 200)
        Rext = R0 + a*np.cos(theta)
        Zext = a*np.sin(theta)
        ax.plot(Rext, Zext, 'k', lw=2)
        
    ax.axis('equal')
    ax.set_xlabel('R [m]', fontsize=14)
    ax.set_ylabel('Z [m]', fontsize=14)
    ax.set_title('Shot #{}@t={}s'.format(shot, time), fontsize=14)
    
    # plot antenna positions and frequencies
    ax.plot([pos[0], pos[0]],[-0.6, 0.6], 'c',
             [pos[0], pos[0]+0.1], [-0.6, -0.6],'c',
             [pos[0], pos[0]+0.1], [0.6, 0.6],'c', lw=2)
    ax.plot([pos[1], pos[1]], [-0.6, 0.6],'y',
             [pos[1], pos[1]+0.1], [-0.6, -0.6],'y',
             [pos[1], pos[1]+0.1], [0.6, 0.6],'y', lw=2)
    ax.plot([pos[2], pos[2]], [-0.6, 0.6],'c--', 
             [pos[2], pos[2]+0.1], [-0.6, -0.6],'c--', 
             [pos[2], pos[2]+0.1], [0.6, 0.6],'c--', lw=2)
             
    ax.text(x=pos[0]+0.1, y=0.5, s='Q1: {} m\n   {} MHz'.format(
            int(pos[0]*100)/100, int(freq[0]*100)/100))
    ax.text(x=pos[1]+0.1, y=0.3, s='Q4: {} m\n   {} MHz'.format(
            int(pos[1]*100)/100, int(freq[1]*100)/100))
    ax.text(x=pos[2]+0.1, y=0.1, s='Q5: {} m\n   {} MHz'.format(
            int(pos[2]*100)/100, int(freq[2]*100)/100))  
    
    # ICRH Resonance theoretical frequencies (H,D)
    R = np.linspace(1.5, 3.5, 501) # radius 30% larger than [R0+/-a]
    B = 0.0073*Itor/R
    B0 = 0.0073*Itor/R0
    fp = e*B/(2*np.pi*m_p*1e6)
    A = 1
    ns = [1,2,3]
    colours = iter(cm.rainbow(np.linspace(0,1,len(ns)*len(freq))))
    

    # Plot the resonance layers for the three harmonics k=1,2,3
    z = np.linspace(-a, a, 101)
    RR, ZZ = np.meshgrid(R, z) 
    for f in freq:
        for n in ns:
            # index at which the resonance condition 
            # starts being satisfied
            ind = np.argwhere(f/n >= fp)
            print(ind)
            if ind.size:
                res_layerH = R[ind[0]]
    
                dummy, R_ripple_max, R_wo_ripple = IC_resonance_radius_ripple(
                                                RR, ZZ, B0, f, n, ep=+1, A=A)    
                dummy, R_ripple_min, R_wo_ripple = IC_resonance_radius_ripple(
                                                RR, ZZ, B0, f, n, ep=-1, A=A)    
      
                print(R_ripple_max)
                
                ax.plot(R_ripple_max, z, 'k')   
                ax.plot(R_ripple_min, z, 'k')   
    
                print('Resonance radius {}H = {} m'.format(n, R_wo_ripple[0]))
                cur_col = next(colours)
                ax.axvline(R_wo_ripple[0], color=cur_col, ls='--', lw=2)
                ax.text(x=res_layerH+0.03, y=0, s='{}H@{} MHz'.format(n,f), 
                        rotation=90, fontsize=14, color=cur_col)
        
if __name__ == '__main__':
    poscouche(freq=55, Itor=900)
    #poscouche(freq=[50, 55, 60], Itor=1250)
    #poscouche(shot=47979, time=10)
       