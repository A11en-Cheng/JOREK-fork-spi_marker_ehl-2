# A small tool used to compare the fluid impurity distribution and the loaded impurity marker
# Usage: python this_scripts jorekxxxxx.h5
# Weikang Tang, weikang.tang@ipp.mpg.de

import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

if __name__ == '__main__':
    # hard-coded param.
    ntor = 15
    # Open the HDF5 file
    file_path = sys.argv[1]
    proj_file = 'projections'+file_path[5:10]+'.h5'
    part_file = 'part'+file_path[5:10]+'.h5'
    with h5py.File(file_path, 'r') as file:
        # List all groups
        print("Keys: %s" % file.keys())
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
    rhon = 0.
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
        Te += np.array(val[7,0,i,:])
        rhon += np.array(val[8,0,i,:])
        rhoimp += np.array(val[9,0,i,:])

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
    
    Te = Te * Tnorm   # Te in eV
    ne = rho + rhoimp*(aux4/10. - 1.)  # ne in JOREK units
    ne = ne * nnorm   # ne in 1e20/m^3
    nn = rhon * nnorm   # nn in 1e20/m^3
    nimp = rhoimp/10.*nnorm  # nimp in 1e20/m^3
    aux5 = aux5*nnorm  # nimp in 1e20/m^3

    # Open the HDF5 file
    with h5py.File(part_file, 'r') as file:
        # Get the data
        coor = file['/groups/001/x'][:]

    RR = np.array(coor[:,0])
    ZZ = np.array(coor[:,1])
    phi = np.array(coor[:,2])
    ind = np.where(RR>0)[0]


    plt.tricontourf(R,Z,nimp,50)
    plt.axis('equal')
    plt.colorbar()
    plt.scatter(RR[ind],ZZ[ind],s=0.5,color='r')
    plt.show()
