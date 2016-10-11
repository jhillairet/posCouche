# -*- coding: utf-8 -*-
import pywed as pw
import numpy as np
    
def formeTS(shot):
    """
    Rparoi,Zparoi,Rext,Zext,t = formeTS(shot)
    
    Arguments:
    - shot: shot number
    
    Returns:
    - R_wall: radius of the vacuum chamber walls [m]
    - Z_wall: height of the vacuum chamber walls [m]
    - R_ext: radius of LCFS [m]
    - Z_ext: height of LCFS [m]
    - t: time [s]
    
    @author: V.Basiuk, J.Hillairet
    """
    if (shot <= 0) or (not isinstance(shot, int)):
        raise ValueError('Shot number should be a positive integer')
        
    # Tore Supra vacuum chamber profile   
    # TODO : WEST vacuum chamber profile
    wall = pw.tsmat(shot, 'APOLO;+Parametres;Paroi')
    R_wall = wall[:,0]
    Z_wall = wall[:,1]
    
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
   
    return R_wall, Z_wall, R_ext, Z_ext, t
    
        
# Below a test code which is run only if this file is executed directly
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    R_wall, Z_wall, R_ext, Z_ext, t = formeTS(47979)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(R_wall, Z_wall, 'k', lw=2)
    ax.axis('equal')
    
    # plasma profile at the middle of the shot
    R_e = R_ext[int(len(R_ext)/2)]
    Z_e = Z_ext[int(len(R_ext)/2)]
    ax.plot(R_e, Z_e, 'b')
    
    