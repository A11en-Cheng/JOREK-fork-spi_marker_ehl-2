#!/usr/bin/python3
import sys
from getopt import *
import numpy as np
import netCDF4 as nc
from scipy.special import legendre

rpos=-1

M=-1
L=-1
nt=1

nz1 = 0
nz2 = 0

outfile=""

try:
    opts, args = getopt(sys.argv[1:],"r:m:l:p:s:o:")
except GetoptError:
    print("Usage: compute_dcoef.py -r <grid point no to integrate at> -m <toroidal modes> -l <poloidal modes> -p <no of toroidal periods> -s <no of grid points to skip> -o <output file> input_file.nc")
    sys.exit(2)

for opt, arg in opts:
    if (opt == "-r"): rpos = int(arg)
    elif (opt == "-m"): M = int(arg)
    elif (opt == "-l"): L = int(arg)
    elif (opt == "-p"): nt = int(arg)
    elif (opt == "-s"): nz1 = int(arg)
    elif (opt == "-o"): outfile = arg

if (M == -1 or L == -1):
    print("Please specify number of toroidal and poloidal modes")
    print("Usage: compute_dcoef.py -r <grid point no to integrate at> -m <toroidal modes> -l <poloidal modes> -p <no of toroidal periods> -s <no of grid points to skip> -o <output file> input_file.nc")
    sys.exit(2)

# 1D integration in the toroidal direction via the trapezoid rule, which
#   reduces to a simple unweighted summation for periodic functions
def integrate_phi(f,dp): return dp*np.sum(f)

# 2D integration via the trapezoid rule in the toroidal direction and 5th order interpolation in the z direction
def integrate2D(f,dz,dp):
    sz = f.shape[0]
    return dz*dp*np.sum(95*(f[0,:] + f[sz-1,:])/288 + 317*(f[1,:] + f[sz-2,:])/240 + 23*(f[2,:] + f[sz-3,:])/30 + 793*(f[3,:] + f[sz-4,:])/720 + 157*(f[4,:] + f[sz-5,:])/160) + dz*dp*np.sum(f[5:(sz-5),:])

def inprodlg(p,q):
    lg = legendre(q)
    zp = np.poly1d(np.zeros(p),r=True)
    intg = (zp*lg).integ()
    return intg(1.0) - intg(-1.0)

ds = nc.Dataset(sys.argv[len(sys.argv)-1])
nr = len(ds.dimensions["rad"])
nz = len(ds.dimensions["zee"])
nf = len(ds.dimensions["phi"])
rmax = float(ds["rmax"][:])
rmin = float(ds["rmin"][:])
zmax = float(ds["zmax"][:])
zmin = float(ds["zmin"][:])
Br = ds["br_001"][:]
Bz = ds["bz_001"][:]
Bp = ds["bp_001"][:]
try:
    insd = ds["inside"][:]
except IndexError:
    insd = np.ones((nf,nz,nr), dtype=np.int)
ds.close()

if (rpos == -1): rpos = int(nr/2)
R0 = rmin + rpos*(rmax - rmin)/(nr - 1)

for i in range(nz):
    if ((insd[:,i,rpos]).all() and (insd[:,nz-i-1,rpos]).all()):
        nz1 += i
        break

Z = zmax - nz1*(zmax - zmin)/(nz - 1)
Z /= R0
nz2 = nz - nz1
nnz = nz2 - nz1

F0 = R0*integrate_phi(Bp[:,int(nz/2),rpos],1.0/nf)

A = np.zeros((L+1,L+1))
b1 = np.zeros(L+1)
b2 = np.zeros(L+1)
b3 = np.zeros(L)
b4 = np.zeros(L)

ac = np.zeros((M+1,L+1))
bc = np.zeros((M+1,L+1))
cc = np.zeros((M+1,L))
dc = np.zeros((M+1,L))

#build A matrix (same for all cases)
for j in range(L+1):
    for k in range(j,L+1):
        A[j,k] = inprodlg(k,j)*Z**k/np.math.factorial(k)

for i in range(M+1):
    sinmp = np.sin(i*nt*np.linspace(0,2*np.pi/nt,nf,endpoint=False))
    cosmp = np.cos(i*nt*np.linspace(0,2*np.pi/nt,nf,endpoint=False))
    for j in range(L+1):
        Plz = legendre(j)(np.linspace(-1.0,1.0,nnz))
        b1[j] = R0*integrate2D((Bp[:,nz1:nz2,rpos]).transpose()*np.outer(Plz,sinmp),2.0/nnz,2.0/nf)/F0
        b2[j] = R0*integrate2D((Bp[:,nz1:nz2,rpos]).transpose()*np.outer(Plz,cosmp),2.0/nnz,2.0/nf)/F0
        if (j < L):
            b3[j] = R0*integrate2D((Br[:,nz1:nz2,rpos]).transpose()*np.outer(Plz,cosmp),2.0/nnz,2.0/nf)/F0
            b4[j] = R0*integrate2D((Br[:,nz1:nz2,rpos]).transpose()*np.outer(Plz,sinmp),2.0/nnz,2.0/nf)/F0
    if (i == 0):
        b3 /= 2
        b4 /= 2
    else:
        ac[i,:] = np.linalg.solve(-i*nt*A, b1)
        bc[i,:] = np.linalg.solve(i*nt*A, b2)
    cc[i,:] = np.linalg.solve(A[:L,:L], b3)
    dc[i,:] = np.linalg.solve(A[:L,:L], b4)

for j in range(L):
    Plz = legendre(j)(np.linspace(-1.0,1.0,nnz))
    b1[j] = R0*integrate2D((Bz[:,nz1:nz2,rpos]).transpose()*np.outer(Plz,np.ones(nf)),2.0/nnz,1.0/nf)/F0

ac[0,1:] = np.linalg.solve(A[:L,:L], b1[:L])

print("R0 = " + str(R0))
print("F0 = " + str(F0))
print(ac)
print(bc)
print(cc)
print(dc)

if (outfile):
    outarr = np.zeros((nt*M+1,L+1,4))
    for i in range(M+1):
        outarr[nt*i,:,0] = ac[i,:]
        outarr[nt*i,:,1] = bc[i,:]
        outarr[nt*i,1:,2] = cc[i,:]
        outarr[nt*i,1:,3] = dc[i,:]
    np.save(outfile, outarr)
