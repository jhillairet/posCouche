# -*- coding: utf-8 -*-
     
import matplotlib.pyplot as plt
import numpy as np

from pywed import tsbase, PyWEDException
from formeTS import formeTS
#from IRFMtb import iround, isempty

from scipy.constants import m_e, m_p, e

def poscouche(shot=None, Itor=None, shot_time=None, freq=55, R0=2.40, a=0.80):
    """
    Poscouche
    
    Calculates and display the locations of Ion Cyclotron Resonance 
    layers for Hydrogen and Deuterium.
    
    Optionnal Arguments:
    - shot: shot number
    - Itor: current in toroidal coils [A]
    - shot_time: default time = time at the middle of the shot  [s]
    - freq: RF frequency [MHz] (Default value 55 MHz)
    - R0: Major radius [m]. (Default value 2.40 m)
    - a: minor radius [m]. (Default value 0.80 m)
    
    Returns:
    - rcouchH
    - rcouchD
    - rco
    - ran
    
    
    @author: Vincent Basiuk, Julien Hillairet
    """      
    # WEST configuration > 48451
    
    # set default values if not passed
    if shot is None:           
        fre = [freq, freq, freq]       
        shot_time = 5
        axe = R0 + 0.1
        pos = None
   
    else: 
        try:
            xa,ta = tsbase(shot, 'samin', nargout=2)
        except PyWEDException:
            xa = 0.8
            ta = 1
            
        try:
            x,t = tsbase(shot, 'srmaj', nargout=2)
        except PyWEDException:
            x = 2.4
            t = 1
            
        try:
            d,t1 = tsbase(shot, 'sd0mag', nargout=2)
        except PyWEDException:
            d = 0.1
            t1 = 1
        
        
        try:
            fre, ant = tsbase(shot, 'SFREQFCI', nargout=2)
        except PyWEDException:
            # TODO : ask user to enter a frequency in MHz                
            #choice = raw_input('sfreqfci signal not available. Please choise a frequency between: ')
            fre = [freq, freq, freq]
            ant = [0, 0, 1]
                   
        # Set the shot_time to the time of max ICRH power
        # Otherwise to the middle of the shot 
        try:
            pfci, t_pfci = tsbase(shot, 'SPUISS', nargout=2)
            shot_time = t_pfci[np.where(pfci == np.max(pfci))]
        except PyWEDException:
            if shot < 28540:
                shot_time = np.max(t)/2
            else:
                shot_time = np.max(ta[xa>0])/2
        
        # Plasma properties at shot_time
        if isinstance(t, np.ndarray):
            R0 = x[np.argmin(np.abs(t - shot_time))]
            a = xa[np.argmin(np.abs(t - shot_time))]
            axe = R0 + d[np.argmin(np.abs(t - shot_time))]
        else:
            R0 = x
            a = xa
            axe = R0 + d
        
        # Antenna positions
        try:
            pos = tsbase(shot, 'SPOSFCI', nargout=1)
        except PyWEDException:
            pos = [3.3, 3.3, 3.3]      
    
    # Current in toroidal coils at the shot_time
    # TODO : may fails for old shot numbers < 20600 ?
    #    if shot < 20100:
    #        Itor = tsmat(shot, 'EXP=T=S;EXP=T=S;ITOR')
    #    if shot < 20600:
    #        Itor = tsmat(shot, 'EXP=T=S;GENERAL;ITOR')
    #    else:
    Itor, tsi = tsbase(shot, 'SITOR', nargout=2)
    Itor = Itor[np.argmin(np.abs(tsi - shot_time))]     
    
    plt.figure()
    # Plot vacuum chamber and plasma LCFS profiles.
    # wall profile available only for shot > 28540        
    if shot > 28540:
        Rparoi, Zparoi, Rext, Zext, t = formeTS(shot)
        ind_t = np.argmin(np.abs(t - shot_time))
        plt.plot(Rparoi, Zparoi, 'k', lw=2)
        plt.plot(Rext[ind_t,:], Zext[ind_t,:], 'b', lw=1)
    else:
        warnings.warn('Old shot number ideal plasma profile plotted')
        theta = np.linspace(0, 2*np.pi, 200)
        Rext = R0 + a*np.cos(theta)
        Zext = a*np.sin(theta)
        plt.plot(Rext, Zext, 'b', lw=1)
        
    plt.axis('equal')
    plt.xlabel('R [m]', fontsize=14)
    plt.ylabel('Z [m]', fontsize=14)
    
    # plot antenna positions and frequencies
    if pos is not None:
        plt.plot([pos[0], pos[0]],[-0.6, 0.6], 'c',
                 [pos[0], pos[0]+0.1], [-0.6, -0.6],'c',
                 [pos[0], pos[0]+0.1], [0.6, 0.6],'c', lw=2)
        plt.plot([pos[1], pos[1]], [-0.6, 0.6],'y',
                 [pos[1], pos[1]+0.1], [-0.6, -0.6],'y',
                 [pos[1], pos[1]+0.1], [0.6, 0.6],'y', lw=2)
        plt.plot([pos[2], pos[2]], [-0.6, 0.6],'c--', 
                 [pos[2], pos[2]+0.1], [-0.6, -0.6],'c--', 
                 [pos[2], pos[2]+0.1], [0.6, 0.6],'c--', lw=2)
                 
        plt.text(x=pos[0]+0.1, y=0.5, s='Q1: {} m'.format(int(pos[0]*100)/100))
        plt.text(x=pos[1]+0.1, y=0.4, s='Q4: {} m'.format(int(pos[1]*100)/100))
        plt.text(x=pos[2]+0.1, y=0.3, s='Q5: {} m'.format(int(pos[2]*100)/100))
        
        plt.text(x=pos[0]+0.5, y=0.5, s=', {} MHz'.format(int(fre[0]*100)/100))
        plt.text(x=pos[1]+0.5, y=0.4, s=', {} MHz'.format(int(fre[1]*100)/100))
        plt.text(x=pos[2]+0.5, y=0.3, s=', {} MHz'.format(int(fre[2]*100)/100))       
    
    # ICRH Resonance theoretical frequencies (H,D)
    # TODO : is m_D really equals to 2 x m_H ?
    B0 = Itor*3.08/1000
    R = np.linspace(R0-a-0.7, R0+a+0.7, 500)
    B = B0*2.37/R
    fp = e*B/(2*np.pi*m_p*1e6)
    fd = e*B/(2*np.pi*m_p*2*1e6)    
    
    # Plot the resonance layers for the three harmonics k=1,2,3
    R_ripple = []    
    zpo = np.linspace(-a, a, 51)
    RR, ZZ = np.meshgrid(R, zpo) 
    for n in [1,2,3]:
        # index at which the resonance condition 
        # starts being satisfied
        ind = np.argwhere(freq/n >= fp)
    
        if ind.size:
            res_layerH = R[ind[0]]
            print('Resonance radius {}H = {} m'.format(n, int(res_layerH*100)/100))
                                
            res_cond, R_ripple, R_wo_ripple = coucherip(RR, ZZ, B0, freq, n, ep=+1, A=1)    
            plt.plot(R_ripple, zpo, 'k')    
            res_cond, R_ripple, R_wo_ripple = coucherip(RR, ZZ, B0, freq, n, ep=-1, A=1)    
            plt.plot(R_ripple, zpo, 'k')   
            print(res_layerH)
            print(R_wo_ripple[0])
            plt.axvline(R_wo_ripple[0], color='b')
            plt.text(x=res_layerH, y=a-0.05, s='{}H'.format(n))
        
if __name__ == '__main__':
    
    poscouche(shot=47979, shot_time=10)
       