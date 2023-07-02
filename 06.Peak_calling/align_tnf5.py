#!/bin/python
import os
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='convert bedpe to tnf5.bed')

parser.add_argument('--in_dir', '-in', type=str, required=True)
parser.add_argument('--out_dir', '-out', type=str, required=True)
parser.add_argument('--replicate', '-rep', type=str, required=True)
parser.add_argument('--celltype', '-celltype', nargs="+", default=None)

args = parser.parse_args()

os.chdir(os.getenv("HOME"))

# List files
in_dir = args.in_dir
out_dir = args.out_dir
replicate = args.replicate

suffix = f'_{replicate}.bedpe'

base_nm = []
in_files = []
for fn in os.listdir(in_dir):
	if fn[-len(suffix):] == suffix:
		base_nm.append(fn[:-len(suffix)])
		in_files.append(f'{in_dir}/{fn}')
		
os.makedirs(out_dir, exist_ok=True)
out_files = [f'{out_dir}/{g}_{replicate}.tnf5.bed' for g in base_nm]

if args.celltype is not None:
	celltype = pd.Index(args.celltype)
	celltype = celltype[celltype.isin(base_nm)]
	print("The following celltypes are in the directory:", celltype.to_list())

	ind = [base_nm.index(g) for g in celltype]
	in_files = np.array(in_files)[ind]
	out_files = np.array(out_files)[ind]

# Run
for inF, outF in zip(in_files, out_files):
	# # read and write
	with open(inF, 'r') as Reader, open(outF, 'w') as Writer:
		for line in Reader:
			lx = line.strip().split("\t")
			lx_a = [lx[i] for i in [0, 1, 2]] + ["N", "1000"] + [lx[8]]
			lx_b = [lx[i] for i in [3, 4, 5]] + ["N", "1000"] + [lx[9]]
			for x in (lx_a, lx_b):
				s = x[5]
				if s == "+":
					x[1] = str(int(x[1]) + 4)
				elif s == "-":
					x[2] = str(int(x[2]) - 5)
				else:
					raise Exception("the string is neither + nor -")

			new_lx = "\t".join(lx_a) + "\n" + "\t".join(lx_b) + "\n"
			Writer.write(new_lx)

	# #
	print(f"\nSaving {outF} successfully")
