# -*- coding: utf-8 -*-
"""
@author: J.Hillairet
"""
import pywed as pw
import numpy as np
    
def vacuum_vessel(shot):
    """
    Get the coordinates of the Tore Supra / WEST vacuum vessel
    
    R_wall, Z_wall = vacuum_vessel(shot)
    
    Arguments:
    - shot: Tore Supra or WEST shot number
    
    Returns:
    - R_wall: radius of the vacuum chamber walls [m]
    - Z_wall: height of the vacuum chamber walls [m]
    
    TODO: once WEST will have started, get the final vacuum vessel coordinates
    """
    if (shot <= 0) or (not isinstance(shot, int)):
        raise ValueError('Shot number should be a positive integer')        
    elif shot < 50000: # Tore Supra vacuum chamber profile      
        wall = pw.tsmat(shot, 'APOLO;+Parametres;Paroi')
        R_wall = wall[:,0]
        Z_wall = wall[:,1]
    else: # WEST vacuum chamber profile
        R_wall, Z_wall = np.loadtxt('WEST_vacuum_vessel.txt', skiprows=1, unpack=True)

    return R_wall, Z_wall
    
    
def LCFS(shot): 
    """
    Get the coordinates of the LCFS as a function of time.
    
    R_ext, Z_ext, t = LCFS(shot)
    
    Arguments:
        shot: Tore Supra or WEST shot number
    
    Returns:
        R_ext: radius of LCFS [m]
        Z_ext: height of LCFS [m]
        t: time [s]
    """
    if (shot <= 0) or (not isinstance(shot, int)):
        raise ValueError('Shot number should be a positive integer')   
        
    # small radius vs time
    y, t = pw.tsbase(shot, 'GRHO', nargout=2)
    t = t[:,0]
    
    # poloidal profile (assumed circular)
    theta = np.arange(0, 24*15, 15) * np.pi/180
    R0 = 2.42
    
    R_ext = R0 + y*np.cos(theta)
    Z_ext = y*np.sin(theta)    
    
    # trick to have a full profile
    R_ext = np.column_stack((R_ext, R_ext[:,0]))
    Z_ext = np.column_stack((Z_ext, Z_ext[:,0]))
       
    return R_ext, Z_ext, t
 
    
# Below a test code which is run only if this file is executed directly
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    R_wall, Z_wall = vacuum_vessel(47979)
    R_ext, Z_ext, t = LCFS(47979)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(R_wall, Z_wall, 'k', lw=2)
    ax.axis('equal')
    
    # plasma profile at the middle of the shot
    R_e = R_ext[int(len(R_ext)/2)]
    Z_e = Z_ext[int(len(R_ext)/2)]
    ax.plot(R_e, Z_e, 'b')
    
    