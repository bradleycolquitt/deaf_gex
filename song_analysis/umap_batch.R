

library(Seurat)
library(hdf5r)
library(tidyverse)
library(qs)
library(future)
library(furrr)

# Directory setup ---------------------------------------------------------


birds = c(
  "pu37bk12",
  "wh32br40",
  "bk27bk53",
  "or74pu64",
  "bu21wh41",
  "bu22wh42",
  "pk90gr56",
  "pk87gr53",
  "or37pu44",
  "or38pu45",
  "gr38gr1",
  "wh21pk76",
  "gr58gr59",
  "gr56gr57",
  "bk88",
  "bk98",
  "bk101",
  "bk89"
  )
data_dirs = map(birds, function(x) sprintf("/home/brad/data2/avgn/%s/clusters", x))

dim_size = 16
sample_size = 1
subsample_size = .5
fname_prefix =  sprintf("syls_sample%s_%s", sample_size, dim_size)

obj_filt_fnames =  map(data_dirs, function(x) file.path(x, fname_prefix, paste0(fname_prefix, "_filtered_seurat.qs")))
out_base_dir = "~/data2/rstudio/birds/singing/song_analysis/avgn"
script_fname = "umap_batch"

# UMAP params -------------------------------------------------------------

dims_to_use = 1:10
umap_name = sprintf("umap%s", max(dims_to_use))

# Run UMAP -------------------------------------------------------------------------

plan(multisession(workers = 9, gc=T))

walk(obj_filt_fnames, function(obj_filt_fname) {
  print(obj_filt_fname)
  obj = qread(obj_filt_fname)

  obj_filt_sub_fname = sub("\\.qs", sprintf("_sample%s.qs", subsample_size), obj_filt_fname)
  
  cells = sample(Cells(obj), size = round(ncol(obj) * subsample_size))
  obj_sub = obj %>%
    subset(cells=cells)
    obj_sub = obj_sub %>%
      FindVariableFeatures() %>%
      RunPCA() %>%
      RunUMAP(dims=dims_to_use, reduction.name=umap_name, verbose = F)
    qsave(obj_sub, obj_filt_sub_fname)
  
  obj_filt_sub_h5_fname = sub("\\.qs", "\\.h5seurat", obj_filt_sub_fname)

  SeuratDisk::SaveH5Seurat(obj_sub, filename=obj_filt_sub_h5_fname, overwrite = TRUE)
})



