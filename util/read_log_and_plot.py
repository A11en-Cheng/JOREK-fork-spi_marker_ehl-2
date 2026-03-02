#!/usr/bin/env python3
""" Search JOREK log file for lines containing "Psi "," u "," zj "," omega "," rho "," T_i "," v_par "," rho_n "," T_e " and 
    plot to evaluate the time evolution of these quantities for numerical diagnostics. The script can be used to extract data for any string, but the default is set to these specific quantities.
"""

import numpy as np
import pylab as plt
import time
import os
import argparse

mydir=os.getcwd()

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-fname', type=str, help="input file name");
parser.add_argument("-text",type=list,default=["Psi "," u "," zj "," omega "," rho "," T_i "," v_par "," rho_n "," T_e "],action="append",help="string to find");
parser.add_argument('-sum',default=False, action='store_true', help="Sum of entries");
parser.add_argument('-n', type=int, help="max number of entries");
parser.add_argument('-s','--save_file', action='store_true', help="save data to file");
 
parser.parse_args(namespace=arg)

if arg.fname!=None:
	fid = os.path.join(mydir,arg.fname)
else:
	fid=os.path.join(mydir,'logfile.out')
	
if arg.text!=None:
	sstr = arg.text
else:
	sstr= ['gmres/solve']

if arg.n!=None:
	nmax = arg.n
else:
	nmax = 0

print(fid)
print(sstr)

dat_min = {s: [] for s in sstr}
dat_max = {s: [] for s in sstr}
with open(fid,encoding = "ISO-8859-1") as infile:
	for line in infile:
		for s in sstr:
			try:
				if s in line:
					parts = line.strip().split()
					if len(parts) >= 2:
						dat_min[s].append(float(parts[-2]))
						dat_max[s].append(float(parts[-1]))
			except Exception as e:
				print("Error reading file: {}".format(e))


counts_min = {s: len(v) for s, v in dat_min.items()}
counts_max = {s: len(v) for s, v in dat_max.items()}
min_len = min([*counts_min.values(), *counts_max.values()]) if (counts_min or counts_max) else 0
print("Number of etries found per key: {}".format(counts_min))
if (nmax==0):
	nmax = min_len
else:
	nmax = min(nmax, min_len)

keys = list(dat_min.keys())
data_min = np.column_stack([np.array(dat_min[k][:nmax]) for k in keys]) if nmax > 0 else np.empty((0, len(keys)))
data_max = np.column_stack([np.array(dat_max[k][:nmax]) for k in keys]) if nmax > 0 else np.empty((0, len(keys)))

# print(data)
if (arg.sum):
	print("Sum min per key = {}".format({k: np.sum(data_min[:, i]) for i, k in enumerate(keys)}))
	print("Sum max per key = {}".format({k: np.sum(data_max[:, i]) for i, k in enumerate(keys)}))
else:
	print("Average min per key = {}".format({k: np.mean(data_min[:, i]) for i, k in enumerate(keys)}))
	print("Average max per key = {}".format({k: np.mean(data_max[:, i]) for i, k in enumerate(keys)}))

if arg.save_file:
	outfile = "logfile_extraction.dat"
	with open(outfile, "w", encoding="utf-8") as out:
		out.write(" ".join([f"{k}_min" for k in keys] + [f"{k}_max" for k in keys]) + "\n")
		if data_min.size and data_max.size:
			merged = np.column_stack([data_min, data_max])
			np.savetxt(out, merged)
	print("Data saved to file: {}".format(outfile))

if data_min.size and data_max.size:
	plt.figure()
	for i, k in enumerate(keys):
		plt.plot(data_min[:, i], label=f"{k} min")
		plt.plot(data_max[:, i], label=f"{k} max")
	plt.ylim(0, 10)
	plt.xlabel("Index")
	plt.ylabel("Value")
	plt.legend()
	plt.tight_layout()
	plt.show()


exit()
