#!/usr/bin/python3
import sys
from getopt import *
import numpy as np
import sympy as sp
import netCDF4 as nc
import field

def integrate(f, Rval, dR, dz, dp, domain):
    return np.sum((f*Rval)[domain.astype(bool)])*dR*dz*dp

R0 = 0.0

M=-1
L=-1
nt=1

outfile=""

try:
    opts, args = getopt(sys.argv[1:],"R:m:l:p:o:")
except GetoptError:
    print("Usage: compute_dcoef.py -R <value of R0> -m <toroidal modes> -l <poloidal modes> -p <no of toroidal periods> -o <output file> input_file.nc")
    sys.exit(2)

for opt, arg in opts:
    if (opt == "-R"): R0 = float(arg)
    elif (opt == "-m"): M = int(arg)
    elif (opt == "-l"): L = int(arg)
    elif (opt == "-p"): nt = int(arg)
    elif (opt == "-o"): outfile = arg

if (M == -1 or L == -1):
    print("Please specify number of toroidal and poloidal modes")
    print("Usage: compute_dcoef.py -R <value of R0> -m <toroidal modes> -l <poloidal modes> -p <no of toroidal periods> -s <no of grid points to skip> -o <output file> input_file.nc")
    sys.exit(2)

domm = field.DommaschkPotentials(np.zeros((1,1,4)))

ds = nc.Dataset(sys.argv[len(sys.argv)-1])
nr = len(ds.dimensions["rad"])
nz = len(ds.dimensions["zee"])
nf = len(ds.dimensions["phi"])
rmax = float(ds["rmax"][:])
rmin = float(ds["rmin"][:])
zmax = float(ds["zmax"][:])
zmin = float(ds["zmin"][:])
BR = ds["br_001"][:]
Bz = ds["bz_001"][:]
Bp = ds["bp_001"][:]
try:
    insd = ds["inside"][:]
except IndexError:
    insd = np.ones((nf,nz,nr), dtype=np.int)
ds.close()

if (R0 == 0.0): R0 = (rmax + rmin)/2

R = np.linspace(rmin, rmax, nr)
z = np.linspace(zmin, zmax, nz)
p = np.linspace(0, 2*np.pi/nt, nf, endpoint=False)
pp, zz, RR = np.meshgrid(p, z, R, indexing="ij")

dim = 1 + 2*M*(L+1) + 2*(M+1)*L
fR = np.zeros((dim, nf,nz,nr))
fz = np.zeros((dim, nf,nz,nr))
fp = np.zeros((dim, nf,nz,nr))

fp[0,:,:,:] = 1.0
i = 1
for l in range(1,L+1):
    fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(domm.D(0,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
    fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(domm.D(0,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
    i += 1
for m in range(1,M+1):
    for l in range(L+1):
        fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
        fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
        fp[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.phi), "numpy")(pp, zz/R0, RR/R0)
        i += 1

for m in range(1,M+1):
    for l in range(L+1):
        fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
        fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
        fp[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.D(m*nt,l).doit(), domm.phi), "numpy")(pp, zz/R0, RR/R0)
        i += 1

for l in range(L):
    fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(domm.N(0,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
    fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(domm.N(0,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
    i += 1
for m in range(1,M+1):
    for l in range(L):
        fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
        fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
        fp[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.cos(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.phi), "numpy")(pp, zz/R0, RR/R0)
        i += 1

for m in range(1,M+1):
    for l in range(L):
        fR[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.R), "numpy")(pp, zz/R0, RR/R0)/R0
        fz[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.Z), "numpy")(pp, zz/R0, RR/R0)/R0
        fp[i,:,:,:] = sp.lambdify((domm.phi, domm.Z, domm.R), sp.diff(sp.sin(m*nt*domm.phi)*domm.N(m*nt,l).doit(), domm.phi), "numpy")(pp, zz/R0, RR/R0)
        i += 1

A = np.zeros((dim,dim))
b = np.zeros(dim)
dR = (rmax - rmin)/(nr-1); dz = (zmax - zmin)/(nz-1); dp = 2*np.pi/(nt*nf)
for i in range(dim):
    b[i] = integrate(BR*fR[i,:,:,:] + Bz*fz[i,:,:,:] + Bp*fp[i,:,:,:]/RR, RR, dR, dz, dp, insd)
    for j in range(dim):
        A[i,j] = integrate(fR[i,:,:,:]*fR[j,:,:,:] + fz[i,:,:,:]*fz[j,:,:,:] + fp[i,:,:,:]*fp[j,:,:,:]/RR**2, RR, dR, dz, dp, insd)

alpha = np.linalg.solve(A, b)
F0 = alpha[0]
ac = np.zeros((M+1,L+1))
bc = np.zeros((M+1,L+1))
cc = np.zeros((M+1,L))
dc = np.zeros((M+1,L))

i = 1
for l in range(1,L+1):
    ac[0,l] = alpha[i]/F0
    i += 1
for m in range(1,M+1):
    for l in range(L+1):
        ac[m,l] = alpha[i]/F0
        i += 1

for m in range(1,M+1):
    for l in range(L+1):
        bc[m,l] = alpha[i]/F0
        i += 1

for m in range(M+1):
    for l in range(L):
        cc[m,l] = alpha[i]/F0
        i += 1

for m in range(1,M+1):
    for l in range(L):
        dc[m,l] = alpha[i]/F0
        i += 1

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
