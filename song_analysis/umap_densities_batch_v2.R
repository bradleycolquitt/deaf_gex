
library(Seurat)
library(hdf5r)
library(tidyverse)
library(qs)
library(future)
library(furrr)
library(cowplot)

source("~/data2/rstudio/birds/utils/db.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
source("~/data2/rstudio/birds/utils/stats.R")


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

obj_filt_fnames =  map(data_dirs, function(x) file.path(x, fname_prefix, paste0(fname_prefix, "_filtered_seurat.qs"))) %>%
  set_names(birds)
out_base_dir = "~/data2/rstudio/birds/singing/song_analysis/avgn"
script_fname = "umap_densities_batch"
out_dir = file.path(out_base_dir, script_fname)
dir.create(out_dir)


# Load metadata -----------------------------------------------------------

info_db = "combined.db"
info_table = "db_info_combined_curated"
info = load_umi_info(info_db, info_table)
info_red = info %>% rename(bird = tags) %>% 
  distinct(procedure2, duration.of.experiment, bird)

# Process data ---------------------------------------------------------------

plan(multisession(workers = 6))

redo = F
umap_to_use = "umap"
if (redo) {
  future_walk(seq_along(obj_filt_fnames), function(i) {
    obj_filt_fname_cur = obj_filt_fnames[[i]]
    bird_cur = birds[i]
    print(bird_cur)
    
    obj = qread(obj_filt_fname_cur)
    
    ### Densities

    umap_dat = Embeddings(obj, umap_to_use) %>% as.data.frame() %>% rownames_to_column()
    group_char = "date_rel"
    group_sym = as.symbol(group_char)
    md = obj@meta.data %>% 
      rownames_to_column()
    umap_dat = umap_dat %>% left_join(md) %>%
      mutate(cluster = !!group_sym)
    
    umap_dat_n = umap_dat %>% group_by(date_rel) %>% summarize(n=n())
    dates = umap_dat_n %>% 
      select(date_rel) %>% arrange() %>%
      unlist() %>% c()
      
    res = map(dates, function(date_to_use) {
      print(date_to_use)
      em_cur = umap_dat %>% filter(date_rel==date_to_use)
      kde = MASS::kde2d(em_cur$UMAP_1, em_cur$UMAP_2, n = 200)
      kde$z
      
    })
    names(res) = dates
    
    res_df = map(res, function(rd) {
      rd_df = as.data.frame(rd) %>%
        rownames_to_column() %>%
        pivot_longer(cols=-rowname, names_to = "axis2", values_to="value") %>%
        mutate(axis2 = as.numeric(sub("V", "", axis2))) %>%
        rename(axis1 = rowname) %>%
        mutate(axis1 = as.numeric(axis1))
      rd_df
    }) %>% set_names(dates) %>% bind_rows(.id="date_rel") %>%
      mutate(date_rel = as.integer(date_rel))
    
    
    baselines = as.character(dates[dates<=0])
    res_baseline_mean = apply(simplify2array(res[baselines]), 1:2, mean)
    res_baseline_sd = apply(simplify2array(res[baselines]), 1:2, sd)
    res_diff = map(res, function(x) (x - res_baseline_mean))
    
    res_diff_df = map(res_diff, function(rd) {
      rd_df = as.data.frame(rd) %>%
        rownames_to_column() %>%
        pivot_longer(cols=-rowname, names_to = "axis2", values_to="value") %>%
        mutate(axis2 = as.numeric(sub("V", "", axis2))) %>%
        rename(axis1 = rowname) %>%
        mutate(axis1 = as.numeric(axis1))
      rd_df
    }) %>% set_names(dates) %>% bind_rows(.id="date_rel") %>%
      mutate(date_rel = as.integer(date_rel))
    
    qsave(res_df, file.path(out_dir, sprintf("%s_%s_density.qs", bird_cur, umap_to_use)))
    qsave(res_diff_df, file.path(out_dir, sprintf("%s_%s_density_diff.qs", bird_cur, umap_to_use)))
  }, .progress=T)
  
  res_dfs = map_dfr(seq_along(obj_filt_fnames), function(i) {
    obj_filt_fname_cur = obj_filt_fnames[[i]]
    bird_cur = birds[i]
    
    res_df = qread(file.path(out_dir, sprintf("%s_%s_density.qs", bird_cur, umap_to_use)))
    res_df$bird = bird_cur
    res_df
  })
  
  res_diff_dfs = map_dfr(seq_along(obj_filt_fnames), function(i) {
    obj_filt_fname_cur = obj_filt_fnames[[i]]
    bird_cur = birds[i]
    
    res_diff_df = qread(file.path(out_dir, sprintf("%s_%s_density_diff.qs", bird_cur, umap_to_use))) 
    res_diff_df$bird = bird_cur
    res_diff_df
  })
} else {
  res_dfs = map_dfr(seq_along(obj_filt_fnames), function(i) {
    obj_filt_fname_cur = obj_filt_fnames[[i]]
    bird_cur = birds[i]
    
    res_df = qread(file.path(out_dir, sprintf("%s_%s_density.qs", bird_cur, umap_to_use)))
    res_df$bird = bird_cur
    res_df
  })
  
  res_diff_dfs = map_dfr(seq_along(obj_filt_fnames), function(i) {
    obj_filt_fname_cur = obj_filt_fnames[[i]]
    bird_cur = birds[i]
    
    res_diff_df = qread(file.path(out_dir, sprintf("%s_%s_density_diff.qs", bird_cur, umap_to_use))) 
    res_diff_df$bird = bird_cur
    res_diff_df
  })
}

# Plotting ----------------------------------------------------------------

walk(birds, function(bird_cur) {
  res_df = res_dfs %>% filter(bird==bird_cur)
  res_diff_df = res_diff_dfs %>% filter(bird==bird_cur)
  dates = unique(res_df$date_rel)
ncols = 4
nrows = ceiling(length(dates) / ncols)
gg = ggplot(res_df, aes(axis1, axis2, fill=value)) +
  geom_raster() +
  facet_wrap(~date_rel, nrow = nrows, ncol= ncols) +
  scale_fill_viridis_c() +
  theme_void()
gg
save_plot(file.path(out_dir, sprintf("%s_%s_density.pdf", bird_cur, umap_to_use)), gg, base_height=2, base_asp=1, ncol=ncols, nrow=nrows)

gg = ggplot(res_diff_df, aes(axis1, axis2, fill=value)) +
  geom_raster() +
  facet_wrap(~date_rel, nrow = nrows, ncol= ncols) +
  scale_fill_gradient2(low=muted("blue"), high=muted("red")) +
  theme_void()
gg
save_plot(file.path(out_dir, sprintf("%s_%s_density_diff.pdf", bird_cur, umap_to_use)), gg, base_height=2, base_asp=1, ncol=ncols, nrow=nrows)

plot(density(res_diff_df$value))

})

# Summary -----------------------------------------------------------------



q99 = function(x) quantile(x, .99)
q01 = function(x) quantile(x, .01)
res_diff_stat = res_diff_dfs %>%
  group_by(bird) %>% 
  filter(value > 0) %>% 
  mutate(value_z = zscore_my(value, baseline=date_rel<=0)) %>%
  group_by(date_rel) %>%
  summarize_at(vars(value, value_z), list(mean=mean, sd=sd, q99=q99, q01=q01))

res_diff_sum = res_diff_dfs %>%
  group_by(bird) %>%
  filter(value > 0) %>% 
  filter(date_rel > -5) %>%
  group_by(date_rel, bird) %>%
  summarize_at(vars(value), list(mean=mean, sd=sd, q99=q99, q01=q01, sum=sum)) %>%
  ungroup() %>%
  mutate(value_z = zscore_my(sum, baseline=date_rel<=0)) %>%
  filter(value_z < 20)

res_diff_sum = res_diff_sum %>%
  left_join(info_red)

x_scale = scale_x_continuous(breaks=seq(min(res_diff_sum$date_rel), max(res_diff_sum$date_rel), 2))

gg = ggplot(res_diff_sum, aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=value_z, group=bird)) + 
  geom_vline(xintercept=0, linetype=2) + 
  theme_bw() + 
  scale_color_manual(values=deaf_colors)
print(gg)

save_plot(file.path(out_dir, sprintf("%s_density_diff_sum_z.pdf", umap_to_use)), gg, base_height=2, base_asp=1.6, ncol=1, nrow=1)


gg = ggplot(res_diff_sum %>% filter(duration.of.experiment==14), aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=value_z, group=bird)) + 
  geom_vline(xintercept=0, linetype=2) + 
  theme_bw() + 
  scale_color_manual(values=deaf_colors)
print(gg)

save_plot(file.path(out_dir, sprintf("%s_density_diff_sum_z_just14.pdf", umap_to_use)), gg, base_height=2, base_asp=1.6, ncol=1, nrow=1)

gg = ggplot(res_diff_sum, aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=value_z, group=bird)) + 
  geom_vline(xintercept=0, linetype=2) + 
  facet_wrap(~duration.of.experiment, ncol=3) +
  theme_bw() + 
  scale_color_manual(values=deaf_colors)
print(gg)

save_plot(file.path(out_dir, sprintf("%s_density_diff_sum_z_duration-split.pdf", umap_to_use)), gg, base_height=2, base_width=3, ncol=3, nrow=1)

gg = ggplot(res_diff_sum, aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=sum, group=bird)) + 
  geom_vline(xintercept=0, linetype=2) + 
  facet_wrap(~duration.of.experiment, ncol=3) +
  theme_bw() + 
  scale_color_manual(values=deaf_colors) + 
  x_scale
print(gg)

save_plot(file.path(out_dir, sprintf("%s_density_diff_sum_duration-split.pdf", umap_to_use)), gg,  base_height=2, base_width=3, ncol=3, nrow=1)


res_diff_sum_1 = res_diff_sum %>%
  group_by(procedure2, date_rel) %>%
  mutate_at(vars(sum, value_z), list(mean=mean, sem=sem, sd=sd)) %>%
  mutate(value_z_mean_up = value_z_mean + value_z_sem,
         value_z_mean_down = value_z_mean - value_z_sem,
         sum_mean_up = sum_mean + sum_sem,
         sum_mean_down = sum_mean - sum_sem)


gg = ggplot(res_diff_sum_1, aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=value_z_mean)) + 
  geom_ribbon(aes(ymin=value_z_mean_down, ymax = value_z_mean_up, fill=procedure2), alpha=.1, color=NA) + 
  geom_vline(xintercept=0, linetype=2) + 
  theme_bw() + 
  scale_color_manual(values=deaf_colors) + 
  scale_fill_manual(values=deaf_colors) + 
  x_scale
print(gg)

save_plot(file.path(out_dir, sprintf("%s_density_sum_z_mean.pdf", umap_to_use)), gg, base_height=2, base_asp=2, ncol=1, nrow=1)

## Fig. 1G
gg = ggplot(res_diff_sum_1, aes(date_rel, color=procedure2)) + 
  geom_line(aes(y=sum_mean)) + 
  geom_ribbon(aes(ymin=sum_mean_down, ymax = sum_mean_up, fill=procedure2), alpha=.1, color=NA) + 
  geom_vline(xintercept=0, linetype=2) + 
  theme_bw() + 
  scale_color_manual(values=deaf_colors) + 
  scale_fill_manual(values=deaf_colors) + 
  x_scale
print(gg)
save_plot(file.path(out_dir, sprintf("%s_density_sum_mean.pdf", umap_to_use)), gg, base_height=2, base_asp=2, ncol=1, nrow=1)

# Difference -- stats ------------------------------------------------------------

library(broom)

diff_stat = res_diff_sum_1 %>% 
    group_by(date_rel) %>%
    do({
      d = .
      print(d$date_rel[1])
      x = d$sum[d$procedure2 == "deaf"]
      y = d$sum[d$procedure2 == "intact"]
      tidy(t.test(x, y))
    })