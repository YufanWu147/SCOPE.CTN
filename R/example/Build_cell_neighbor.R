library(tidyverse)

fig_path = "figures/" # change this line to figure directory

cell_neighbor_results_maxdist = pbmcapply::pbmclapply(
  1:length(slide_list), function(i) {
    Build_cell_neighbor_maxdist(
      dat_slide=dat_in[which(dat_in$ImageID==slide_list[i]), ],
      knn_num = 20, cell_types_available = cell_types_available,
      min.radius = 100, max.radius = 500, min.n.neighbors = 10)},
  mc.cores = 1)

# Visualize the distribution of the 99th percentile of the distances for
# selecting knn_num, min.radius and max.radius
dists_99th_perc <- lapply(cell_neighbor_results_maxdist, "[[", "dist.99th.perc") %>%
  `names<-`(slide_list) %>%
  as.data.frame %>% t()

hist(dists_99th_perc)

# Visualize cells with less than min.n.neighbors cells in their neighborhoods
cell_neighbor_table_maxdist = lapply(
  cell_neighbor_results_maxdist, "[[", "proportions_df") %>%
  do.call("rbind", .)

na_rows = apply(cell_neighbor_table_maxdist, 1, function(x) {any(is.na(x))})
dat_in %>%
  mutate(is_NaN = ifelse(CellID %in% names(na_rows[na_rows]), "Yes", "No")) %>%
  arrange(is_NaN) %>%
  ggplot(aes(x = X, y = Y, color = is_NaN)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey90", "red")) +
  facet_wrap( ~ ImageID) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void(base_size = 10) +
  coord_fixed(expand = FALSE)
ggsave(paste0(fig_path, "max_dist_cells_less_10_nns.png"), width = 10, height = 7)

