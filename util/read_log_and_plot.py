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

dat1 = {s: [] for s in sstr}
with open(fid,encoding = "ISO-8859-1") as infile:
	for line in infile:
		for s in sstr:
			try:
				if s in line:
					line = line.strip().split()
					dat1[s].append(float(line[-1]))
			except Exception as e:
				print("Error reading file: {}".format(e))


counts = {s: len(v) for s, v in dat1.items()}
min_len = min(counts.values()) if counts else 0
print("Number of etries found per key: {}".format(counts))
if (nmax==0):
	nmax = min_len
else:
	nmax = min(nmax, min_len)

keys = list(dat1.keys())
data = np.column_stack([np.array(dat1[k][:nmax]) for k in keys]) if nmax > 0 else np.empty((0, len(keys)))

# print(data)
if (arg.sum):
	print("Sum per key = {}".format({k: np.sum(data[:, i]) for i, k in enumerate(keys)}))
else:
	print("Average per key = {}".format({k: np.mean(data[:, i]) for i, k in enumerate(keys)}))

if arg.save_file:
	outfile = "logfile_extraction.dat"
	with open(outfile, "w", encoding="utf-8") as out:
		out.write(" ".join(keys) + "\n")
		if data.size:
			np.savetxt(out, data)
	print("Data saved to file: {}".format(outfile))

if data.size:
	plt.figure()
	for i, k in enumerate(keys):
		plt.plot(data[:, i], label=k)
	plt.ylim(0, 10)
	plt.xlabel("Index")
	plt.ylabel("Value")
	plt.legend()
	plt.tight_layout()
	plt.show()


exit()
