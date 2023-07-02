import os
import argparse
import gzip
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='split bedpe based on celltype')

parser.add_argument('--donor', '-d', type=str, required=True)
parser.add_argument('--batch', '-b', type=str, required=True)

parser.add_argument('--input_file', '-in', type=str, required=True)
parser.add_argument('--annot_file', '-ano', type=str, required=True)
parser.add_argument('--out_dir', '-out', type=str, required=True)

parser.add_argument('--celltype', '-celltype', nargs="+", default=None)

args = parser.parse_args()

os.chdir(os.getenv("HOME"))

# Load Data
cell_annot = pd.read_csv(args.annot_file)
if args.celltype is not None:
    cell_annot = cell_annot[cell_annot["celltype"].isin(args.celltype)]

print("Split for the following celltypes:", cell_annot["celltype"].unique())	

celltype_df = pd.Series(
    cell_annot["celltype"].to_numpy(), index=cell_annot["uniq_barcode"]
)
celltype_dict = celltype_df[celltype_df.index.sort_values()].to_dict()

# output file
os.makedirs(args.out_dir, exist_ok=True)

outf_dict = dict()
for g in cell_annot["celltype"].unique():
    outf_dict[g] = f'{args.out_dir}/{g}_{args.donor}.bedpe'

# Run
# # gzip open in text mode: rt
bedpeF = gzip.open(args.input_file, 'rt')

write_dict = dict()
for k in outf_dict:
    write_dict[k] = open(outf_dict[k], 'w')

# # write
new_barcode_list = []
avail_clu = []
previous_old = ""

for line in bedpeF:
    lx = line.strip().split("\t")
    read_name = lx[6]
    old_barcode = read_name.split(":")[0]

    if old_barcode != previous_old:
        previous_old = old_barcode
        new_barcode = "_".join((args.batch, old_barcode))
        exist = new_barcode in celltype_dict
        if exist:
            new_barcode_list.append(new_barcode)
            clu = celltype_dict[new_barcode]
            avail_clu.append(clu)
    # inherit the previous "exist"
    elif exist:
        lx[6] = new_barcode
        wline = "\t".join(lx) + "\n"
        write_dict[clu].write(wline)

# # close
bedpeF.close()
for k in write_dict:
    write_dict[k].close()

# # delete empty file
avail_clu = np.unique(avail_clu)
all_clu = pd.Index(outf_dict.keys())
empty_clu = all_clu[~all_clu.isin(avail_clu)]

for k in empty_clu:
    os.remove(outf_dict[k])

# # cell numbers
ncells = pd.Series(new_barcode_list).nunique()
print(f'Split successfully. There are {ncells} cells.')
