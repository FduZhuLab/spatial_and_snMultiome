# Introduction

The scripts for analysis of **"Gene regulatory innovations underlying human neocortical diversification"**.

## 01. Clustering

Clustering analysis including preprocessing and subtype identification by Seurat.

## 02. DEG

Identification of area-specific DEGs at the subclass level.

## 03. Integration

Integration with the public reference datasets to determine the cell-type or compare with other species.

## 04. Spatial Deconvolution

Deconvolution of the spatial transcriptome data, including nucleus segmentation and the subclass/subtype weight prediction.

## 05. CellChat

Cell-cell interaction analysis with cellchat.

## 06. Peak Calling

Generation of reproducible ATAC peaks according to the pipeline of [Bing Ren et al., 2021](https://www.nature.com/articles/s41586-021-03604-1).

## 07. Denosing

Denoising of ATAC data using cisTopic.

## 08. DAP

Identification of area-specific differential accessible peaks (DAPs) at the subclass level.

## 09. AUP Enrichment

Enrichment analysis of AUPs in evolutionary ages and TE families using Fisher's exact test.

## Map files

These are the files that map the PBS_ARRAYID to a specific batch during job submission.
