# append batch+barcode
1. `append_CB.sh`: atac_possorted_bam.bam --> possorted.CB.bam

```bash
qsub append_CB.sh -t 1-8
```

# deduplicated and name-sorted
2. `filter_bam.sh`: possorted.CB.bam --> filtered_dedup_nsrt.bam

```bash
qsub filter_bam.sh -t 1-8
```

# convert to bedpe
3. `bam2bedpe.sh`: filtered_dedup_nsrt.bam --> bedpe.gz

```bash
qsub bam2bedpe.sh -t 1-8
```

# split bedpe by group, and generate pseudo replicate

the group is subclasses in batch (each donor replicate of each area)

4. `split_bedpe.sh` + `split_bedpe.py`: bedpe.gz --> *area*\_*subclass*\_*replicate*.bedpe

	replicate: No_61, No_53

```bash
qsub split_bedpe.sh -t 1-8
```

5. `gen_pseudo.sh`: *area*\_*subclass*\_*replicate*.bedpe --> *area*\_*subclass*\_*pseudorep*.bedpe

	pseudorep: pseudo00, pseudo01

```bash
qsub gen_pseudo.sh -t 1-4
```

# bedpe to bed and correct Tnf5 insertion
6. `align_tnf5.sh` + `align_tnf5.py`: *area*\_*subclass*\_*rep*.bedpe --> *area*\_*subclass*\_*rep*.tnf5.bed

	rep: No_61, No_53, pseudo00, pseudo01

```bash
qsub align_tnf5.sh -t 1-16
```

7. `pool_replicate.sh`: *area*\_*subclass*\_*rep*.tnf5.bed --> *area*\_*subclass*.tnf5.bed

```bash
qsub pool_replicate.sh -t 1-4
```

# peak calling
8. `call_macs2.sh`: *area*\_*subclass*\_*rep*.tnf5.bed --> macs2 results; `call_macs2_pool.sh`: *area*\_*subclass*.tnf5.bed --> macs2 results

```bash
qsub call_macs2.sh -t 1-16 && qsub call_macs2_pool.sh -t 1-4
```

# reproducible peaks
9. `overlap.sh`: macs2 results of pool, replicate and pseudorep --> overlapped peak and its summit 

	**remove your old summit list and old peak list if they are exists**

```bash
qsub overlap.sh -t 1-4
```

# Merge all
10. `merge_peak.sh` + `merge_peak.R`: list of area-subclasses' peak summit --> final.union.bed

```bash
qsub merge_peak.sh
```

# modify fragments file with new barcode
`frag_rename.sh`: atac_fragments.tsv.gz --> atac_fragments.tsv.gz (with new barocde)

```bash
qsub frag_rename.sh -t 1-8
```

# count peaks with *Signac*

using Signac (R package): `count_atac.sh` + `count_atac.R`
