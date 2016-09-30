# -*- coding: utf-8 -*-
import pywed as pw
import numpy as np
    
def formeTS(pulse):
    """
    Rparoi,Zparoi,Rext,Zext,t = formeTS(numchoc)
    
    Arguments
    - pulse: pulse number
    
    Returns:
    - R_wall: radius of the vacuum chamber walls [m]
    - Z_wall: height of the vacuum chamber walls [m]
    - R_ext: radius of LCFS [m]
    - Z_ext: height of LCFS [m]
    - t: time [s]
    
    @author: V.Basiuk, J.Hillairet
    """
    # vacuum chamber profile    
    wall = pw.tsmat(pulse, 'APOLO;+Parametres;Paroi')
    R_wall = wall[:,0]
    Z_wall = wall[:,1]
    
    # small radius vs time
    y, t = pw.tsbase(pulse, 'GRHO', nargout=2)
    t = t[:,0]
    
    # poloidal profile (assumed circular)
    theta = np.arange(0, 24*15, 15) * np.pi/180
    R0 = 2.42
    
    R_ext = R0 + y*np.cos(theta)
    Z_ext = y*np.sin(theta)    
    
    # trick to have a full profile
    R_ext = np.column_stack((R_ext, R_ext[:,0]))
    Z_ext = np.column_stack((Z_ext, Z_ext[:,0]))
   
    return R_wall, Z_wall, R_ext, Z_ext, t
    

# Below a test code which is run only if this file is executed directly
if __name__ == '__main__':
    from matplotlib.pyplot import *
    
    R_wall, Z_wall, R_ext, Z_ext, t = formeTS(47979)
    plot(R_wall, Z_wall, 'k', lw=2)
    axis('equal')
    
    # plasma profile at the middle of the pulse
    R_e = R_ext[np.round(R_ext.shape[0]/2)]
    Z_e = Z_ext[np.round(R_ext.shape[0]/2)]
    plot(R_e, Z_e, 'b')
    
    