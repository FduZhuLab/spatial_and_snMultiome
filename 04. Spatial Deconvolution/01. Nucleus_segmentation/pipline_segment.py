from stardist.models import StarDist2D
from csbdeep.data import Normalizer, normalize_mi_ma

import os
import gc
import argparse
import functools

import numpy as np
import pandas as pd

import scanpy as sc
import squidpy as sq

import tensorflow as tf

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

os.chdir(os.getenv("HOME"))

print = functools.partial(print, flush=True)

# Argements
num_cores = int(os.getenv("SLURM_CPUS_PER_TASK"))

parser = argparse.ArgumentParser()
parser.add_argument('--donor', '-d', type=str, required=True)
parser.add_argument('--region', '-r', type=str, required=True)
parser.add_argument('--spt_base', '-spt', type=str, required=True)
parser.add_argument('--img_base', '-img', type=str, required=True)
parser.add_argument('--output_base', '-o', type=str, required=True)
parser.add_argument('--suffix', '-suffix', nargs='?', type=str, default="")

parser.add_argument('--prob_thresh', nargs='?', type=float, default=0.692478)
parser.add_argument('--nms_thresh', nargs='?', type=float, default=0.3)
parser.add_argument('--block_size', nargs='?', type=int, default=3000)

# parser.add_argument('--n_jobs', nargs='?', type=int, default=num_cores)

args = parser.parse_args()

# Load Data
sample_dir = '/'.join([args.spt_base, args.donor, args.region])
spt_ada = sc.read('/'.join([sample_dir, "spt_adata.h5ad"]))

spt_img = sq.im.ImageContainer()
spt_img.add_img(
    '/'.join([args.img_base, args.donor, args.region, "image.jpg"]),
    layer="image"
)

# Test GPU
print("\ntensorflow is bulit with cuda:", tf.test.is_built_with_cuda())
print("using devices:", tf.config.experimental.list_physical_devices("GPU")[0])
print("GPU memory usage:", tf.config.experimental.get_memory_info("GPU:0"))

# cut blank area
def cut_bdry(img, adata, key="spatial"):
    position = adata.obsm[key]
    key_ = list(adata.uns[key].keys())[0]
    spot_diam = adata.uns[key][key_]["scalefactors"]["spot_diameter_fullres"]
    
    x_u = position[:, 0].max()
    x_l = position[:, 0].min()
    y_u = position[:, 1].max()
    y_l = position[:, 1].min()
    
    center_x = int((x_u + x_l)/2)
    center_y = int((y_u + y_l)/2)
    
    size_x = int((x_u - x_l)/2 + spot_diam)
    size_y = int((y_u - y_l)/2 + spot_diam)

    crop = img.crop_center(y=center_y, x=center_x, radius=(size_y, size_x))
    return crop


spt_img = cut_bdry(spt_img, spt_ada)
# spt_ada = spt_img.subset(spt_ada)

# define stardist
class MyNormalizer(Normalizer):
    def __init__(self, mi, ma):
            self.mi, self.ma = mi, ma
    def before(self, x, axes):
        return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)
    def after(*args, **kwargs):
        assert False
    @property
    def do_after(self):
        return False
    

def stardist_2D_he_large(
    img, block_size=4096, min_overlap=128, context=None, n_tiles=None,
    normalizer=None, fold=8, nms_thresh=None, prob_thresh=None
):
    # The crucial assumption is that all predicted object instances are smaller than
    # the provided `min_overlap`. Also, it must hold that: min_overlap + 2*context < block_size    
    if context is None:
        context = min_overlap
    if min_overlap + 2*context >= block_size:
        raise Exception('block_size is too small')
 
    # Make sure to normalize the input image beforehand
    # or supply a normalizer to the prediction function.
    if normalizer is None:
        mi, ma = np.percentile(img[::fold], [1, 99.8])
        normalizer = MyNormalizer(mi, ma)

    model = StarDist2D.from_pretrained('2D_versatile_he')
    labels, _ = model.predict_instances_big(
        img, normalizer=normalizer, axes='YXC',
        block_size=block_size, min_overlap=min_overlap,
        context=context, n_tiles=n_tiles,
        nms_thresh=nms_thresh, prob_thresh=prob_thresh
    )
    return labels


# Run Stardist
print(
    "\n", "block_size is:", args.block_size, "prob_thresh is:", args.prob_thresh,
    "nms_thresh is:", args.nms_thresh, "\n"
)

gc.collect()
sq.im.segment(
    img=spt_img, layer="image", channel=None,
    method=stardist_2D_he_large,
    layer_added='segmented_stardist',
    block_size=args.block_size,
    prob_thresh=args.prob_thresh,
    nms_thresh=args.nms_thresh
)

print(f"\nNumber of segments in crop: {len(np.unique(spt_img['segmented_stardist']))}")

gc.collect()

# save
img_out_dir = '/'.join([args.output_base, args.donor, args.region, "zarr" + args.suffix])
os.makedirs(img_out_dir, exist_ok=True)
spt_img.save(img_out_dir)

print("\nSegmentation Success")
