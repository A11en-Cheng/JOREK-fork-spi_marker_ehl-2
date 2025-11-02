# A small tool to plot the fields directly from the JOREK restart file
# Plots are based on the simulation grid, i.e. not using n_sub to increase the resolution
# Usage: python this_scripts jorekxxxxx.h5
# Weikang Tang, weikang.tang@ipp.mpg.de

import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

#plt.style.use('seaborn-v0_8')
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.weight'] = 'bold'
font = {'family': 'Serif',
'style': 'normal',
'weight': 'bold',
'color': 'black',
'size': 13,
}

if __name__ == '__main__':
    # hard-coded param.
    ntor = 15
    # Open the HDF5 file
    file_path = sys.argv[1]
    proj_file = 'projections'+file_path[5:10]+'.h5'
    with h5py.File(file_path, 'r') as file:
        # List all groups
#        print("Keys: %s" % file.keys())
        # Get the data
        coor = file['x'][:]
        val = file['values'][:]

    psi = 0.
    U = 0.
    zj = 0.
    w = 0.
    rho = 0.
    Ti = 0.
    v_par = 0.
    rhoimp = 0.
    Te = 0.
    aux1 = 0.
    aux2 = 0.
    aux3 = 0.
    aux4 = 0.
    aux5 = 0.
    R = np.array(coor[0,0,0,:])
    Z = np.array(coor[1,0,0,:])
    for i in range(ntor):
        psi += np.array(val[0,0,i,:])
        U += np.array(val[1,0,i,:])
        zj += np.array(val[2,0,i,:])
        w += np.array(val[3,0,i,:])
        rho += np.array(val[4,0,i,:])
        Ti += np.array(val[5,0,i,:])
        v_par += np.array(val[6,0,i,:])
        rhoimp += np.array(val[7,0,i,:])
        Te += np.array(val[8,0,i,:])

    # open projection file
    try:
        with h5py.File(proj_file,'r') as file:
            proj = file['values'][:]
        for i in range(ntor):
            aux1 += np.array(proj[0,0,i,:])
            aux2 += np.array(proj[1,0,i,:])
            aux3 += np.array(proj[2,0,i,:])
            aux4 += np.array(proj[3,0,i,:])
            aux5 += np.array(proj[4,0,i,:])
    except:
        print('warning: no projection file imported')

    #some normalization factors
    qe = 1.6e-19
    mu0 = 4e-7*np.pi
    nnorm = 0.885173e20
    Tnorm = 1./qe/mu0/nnorm
    rhonorm = nnorm*2.*1.673e-27
    tnorm = 6.1e-7
    jnorm = 1e-6/mu0
    
    Te = Te * Tnorm   # Te in eV
    ne = rho + rhoimp*(aux4/10. - 1.)  # ne in JOREK units
    ne = ne * nnorm   # ne in 1e20/m^3
    nimp = rhoimp/10. * nnorm  # nimp in 1e20/m^3 
    v_par = v_par/tnorm*1e-3*1.8    # vpar in km/s  ! Btot is approximated by Bt 
    zj = zj*jnorm    # zj in MA/m^2
    U = U/tnorm   # U in m/s

    fig, ax = plt.subplots(2, 3, figsize=(20, 10))
    tcf0 = ax[0,0].tricontourf(R,Z,Te,50,cmap='hot')
    ax[0,0].axis('equal')
    ax[0,0].set_ylabel('Z [m]',fontdict=font)
    ax[0,0].set_xlabel('R [m]',fontdict=font)
    ax[0,0].set_title('Te [eV]',fontdict=font)
    fig.colorbar(tcf0)
    tcf1 = ax[0,1].tricontourf(R,Z,nimp,50)
    ax[0,1].axis('equal')
    ax[0,1].set_ylabel('Z [m]',fontdict=font)
    ax[0,1].set_xlabel('R [m]',fontdict=font)
    ax[0,1].set_title(r'$nimp\ [m^{-3}]$',fontdict=font)
    fig.colorbar(tcf1)
    tcf2 = ax[0,2].tricontourf(R,Z,v_par,50,cmap='plasma')
    ax[0,2].axis('equal')
    ax[0,2].set_ylabel('Z [m]',fontdict=font)
    ax[0,2].set_xlabel('R [m]',fontdict=font)
    ax[0,2].set_title('vpar [km/s]',fontdict=font)
    fig.colorbar(tcf2)
    tcf3 = ax[1,0].tricontourf(R,Z,zj,50,cmap='PiYG')
    ax[1,0].axis('equal')
    ax[1,0].set_ylabel('Z [m]',fontdict=font)
    ax[1,0].set_xlabel('R [m]',fontdict=font)
    ax[1,0].set_title(r'$zj\ [MA/m^2]$',fontdict=font)
    fig.colorbar(tcf3)
    tcf4 = ax[1,1].tricontourf(R,Z,U,50,cmap='seismic')
    ax[1,1].axis('equal')
    ax[1,1].set_ylabel('Z [m]',fontdict=font)
    ax[1,1].set_xlabel('R [m]',fontdict=font)
    ax[1,1].set_title('U [m/s]',fontdict=font)
    fig.colorbar(tcf4)
    tcf5 = ax[1,2].tricontourf(R,Z,w,50,cmap='jet')
    ax[1,2].axis('equal')
    ax[1,2].set_ylabel('Z [m]',fontdict=font)
    ax[1,2].set_xlabel('R [m]',fontdict=font)
    ax[1,2].set_title('vorticity',fontdict=font)
    fig.colorbar(tcf5)

    plt.tight_layout()
    plt.show()
