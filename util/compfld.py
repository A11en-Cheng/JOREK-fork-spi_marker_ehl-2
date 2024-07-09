#!/usr/bin/python3
############################################################################################################
# This script depends on the Python implementation of Dommaschk potentials from the BOUT++/Zoidberg project,
# which is available at https://github.com/boutproject/zoidberg
# Only the files zoidberg/field.py, zoidberg/boundary.py and zoidberg/progress.py are needed.
# In field.py, the import statements for boundary and update_progress should be replaced with the following:
#    import boundary
#    from progress import update_progress
############################################################################################################
import sys
from getopt import *
import numpy as np
import netCDF4 as nc
import field

nt=1
R0=1.0
F0=1.0
sz=1
nz1=0

outfile=""

try:
    opts, args = getopt(sys.argv[1:],"R:F:p:z:s:o:")
except GetoptError:
    print("Usage: compfld.py -R <value of R0> -F <value of F0> -p <number of toroidal periods> -z <step size in z direction> -s <distance from boundary to ignore (dz)> -o <output file name> coefficients_file.npy vacuum_field_file.nc")
    sys.exit(2)

for opt, arg in opts:
    if (opt == "-R"): R0 = float(arg)
    elif (opt == "-F"): F0 = float(arg)
    elif (opt == "-p"): nt = int(arg)
    elif (opt == "-z"): sz = int(arg)
    elif (opt == "-s"): nz1 = int(arg)
    elif (opt == "-o"): outfile = arg

coef = np.load(sys.argv[len(sys.argv)-2])
fld = field.DommaschkPotentials(coef)

ds = nc.Dataset(sys.argv[len(sys.argv)-1])
nr = len(ds.dimensions["rad"])
nz = len(ds.dimensions["zee"])
nf = len(ds.dimensions["phi"])
rmax = float(ds["rmax"][:])
rmin = float(ds["rmin"][:])
zmax = float(ds["zmax"][:])
zmin = float(ds["zmin"][:])
Br1 = ds["br_001"][:]
Bz1 = ds["bz_001"][:]
Bp1 = ds["bp_001"][:]
try:
    insd = ds["inside"][:]
except IndexError:
    insd = np.ones((nf,nz,nr), dtype=np.int)
ds.close()

R = np.linspace(rmin, rmax, nr)
Z = np.linspace(zmin, zmax, nz)
P = np.linspace(0, 2*np.pi/nt, nf, endpoint=False)

RR, PP, ZZ = np.meshgrid(R, P, Z[::sz], indexing="ij")
Br2 = F0*fld.Bxf(RR/R0,PP,ZZ/R0)/R0
Bz2 = F0*fld.Bzf(RR/R0,PP,ZZ/R0)/R0
Bp2 = F0*fld.Byf(RR/R0,PP,ZZ/R0)/R0

Br2 = np.moveaxis(Br2, 0, -1)
Bz2 = np.moveaxis(Bz2, 0, -1)
Bp2 = np.moveaxis(Bp2, 0, -1)
RR  = np.moveaxis(RR,  0, -1)

intdiff = np.sum((((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)*RR)[insd[:,::sz,:].astype(bool)]*(rmax-rmin)*(zmax-zmin)*2*sz*np.pi/(nr*nz*nt*nf))
intdiff1 = np.sum((((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)*RR/(Br1[:,::sz,:]**2 + Bz1[:,::sz,:]**2 + Bp1[:,::sz,:]**2))[insd[:,::sz,:].astype(bool)]*(rmax-rmin)*(zmax-zmin)*2*sz*np.pi/(nr*nz*nt*nf))
intdiff2 = np.sum((((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)*RR/(Br2**2 + Bz2**2 + Bp2**2))[insd[:,::sz,:].astype(bool)]*(rmax-rmin)*(zmax-zmin)*2*sz*np.pi/(nr*nz*nt*nf))
vol = np.sum(RR[insd[:,::sz,:].astype(bool)]*(rmax-rmin)*(zmax-zmin)*2*sz*np.pi/(nr*nz*nt*nf))
print("<|B2 - B1|^2> = " + str(intdiff/vol))
print("<(|B2 - B1|/|B1|)^2> = " + str(intdiff1/vol))
print("<(|B2 - B1|/|B2|)^2> = " + str(intdiff2/vol))

maxdiff = np.max(((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)[insd[:,::sz,:].astype(bool)])
maxdiff1 = np.max((((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)/(Br1[:,::sz,:]**2 + Bz1[:,::sz,:]**2 + Bp1[:,::sz,:]**2))[insd[:,::sz,:].astype(bool)])
maxdiff2 = np.max((((Br2 - Br1[:,::sz,:])**2 + (Bz2 - Bz1[:,::sz,:])**2 + (Bp2 - Bp1[:,::sz,:])**2)/(Br2**2 + Bz2**2 + Bp2**2))[insd[:,::sz,:].astype(bool)])
print("max(B2 - B1)^2 = " + str(maxdiff))
print("max(|B2 - B1|/|B1|)^2 = " + str(maxdiff1))
print("max(|B2 - B1|/|B2|)^2 = " + str(maxdiff2))

if (outfile):
    ds = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
    ds.createDimension("rad",nr)
    ds.createDimension("zee",len(Z[::sz]))
    ds.createDimension("phi",nf)
    ds.createVariable("br_orig", "f8", ("phi","zee","rad"))
    ds.createVariable("bz_orig", "f8", ("phi","zee","rad"))
    ds.createVariable("bp_orig", "f8", ("phi","zee","rad"))
    ds.createVariable("br_domm", "f8", ("phi","zee","rad"))
    ds.createVariable("bz_domm", "f8", ("phi","zee","rad"))
    ds.createVariable("bp_domm", "f8", ("phi","zee","rad"))
    ds.createVariable("rmax", "f8")
    ds.createVariable("rmin", "f8")
    ds.createVariable("zmax", "f8")
    ds.createVariable("zmin", "f8")
    ds["rmax"][:] = rmax
    ds["rmin"][:] = rmin
    ds["zmax"][:] = zmax
    ds["zmin"][:] = zmin
    ds["br_orig"][:] = Br1[:,::sz,:]
    ds["bz_orig"][:] = Bz1[:,::sz,:]
    ds["bp_orig"][:] = Bp1[:,::sz,:]
    ds["br_domm"][:] = Br2
    ds["bz_domm"][:] = Bz2
    ds["bp_domm"][:] = Bp2
    ds.close()
