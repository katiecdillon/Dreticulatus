setwd("/scratch/kcd88651/ticks/Dermacentor_Reticulatus/pixy/DretUK_merge/pixy_input")

library(tidyverse)
library(ggh4x)
library(ggpubr)

# -----------------------------
# LOAD DATA
# -----------------------------
pixy_to_long <- function(pixy_files) {
  
  single_pop_stats <- c("pi", "watterson_theta", "tajima_d")
  pixy_df <- list()
  
  for (i in seq_along(pixy_files)) {
    
    stat_file_type <- sub(".*_(pi|tajima_d|watterson_theta)\\.txt$",
                          "\\1",
                          basename(pixy_files[i]))
    
    df <- read_delim(pixy_files[i], delim = "\t")
    
    if (stat_file_type %in% single_pop_stats) {
      # In the single_pop_stats branch, remove the as.character() mutate entirely:
      df <- df %>%
        select(-any_of("tajima_d_s_counts")) %>%   # drop non-numeric before pivoting
        pivot_longer(cols = -c(pop, window_pos_1, window_pos_2, chromosome),
                     names_to = "statistic",
                     values_to = "value") %>%
        mutate(value = as.numeric(value)) %>%       # ensure numeric after pivot
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
    } else {
      df <- df %>%
        pivot_longer(cols = -c(pop1, pop2, window_pos_1, window_pos_2, chromosome),
                     names_to = "statistic",
                     values_to = "value") %>%
        mutate(value = as.numeric(value))
    }
    
    pixy_df[[i]] <- df
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_folder <- "/scratch/kcd88651/ticks/Dermacentor_Reticulatus/pixy/DretUK_merge/pixy_input"
pixy_files <- list.files(pixy_folder, pattern = "\\.txt$", full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

# -----------------------------
# RENAME CHROMOSOMES
# -----------------------------
chromosome_key <- data.frame(chr_numb = c(1,2,3,4,5,6,7,8,9,10,11),
                             chr_name = c("CM170536.1","CM170537.1",
                                          "CM170538.1","CM170539.1",
                                          "CM170540.1","CM170541.1",
                                          "CM170542.1","CM170543.1",
                                          "CM170544.1","CM170545.1",
                                          "CM170546.1"))

pixy_df <- pixy_df %>%
  left_join(chromosome_key, by = c("chromosome" = "chr_name")) %>%
  mutate(chromosome = chr_numb) %>%
  select(-chr_numb)

write_tsv(pixy_df, "DretUKmerge_pixy_df.tsv")

sink("DretUKmerge_pixy_checks.log")

cat("Check 1: NAs per statistic\n")
pixy_df %>%
  filter(statistic %in% c("avg_pi","avg_watterson_theta","avg_dxy","avg_hudson_fst","tajima_d")) %>%
  group_by(statistic) %>%
  summarise(n_na = sum(is.na(value))) %>%
  print()

cat("\nCheck 2: Value limits per statistic\n")
pixy_df %>%
  filter(statistic %in% c("avg_pi","avg_watterson_theta","avg_dxy","avg_hudson_fst","tajima_d")) %>%
  group_by(statistic) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ) %>%
  print()

sink()

# -----------------------------
# PLOT
# -----------------------------
pixy_labeller <- as_labeller(c(avg_pi              = "pi",
                               avg_dxy             = "D[XY]",
                               avg_wc_fst          = "F[ST]",
                               avg_hudson_fst      = "F[ST]",
                               avg_watterson_theta = "theta[W]",
                               tajima_d            = "Tajima*minute*s~D"),
                             default = label_parsed)

# -----------------------------
# PI
# -----------------------------
p1 <- pixy_df %>%
  mutate(chrom_color_group = if_else(chromosome %% 2 == 1, "even", "odd")) %>%
  filter(statistic == "avg_pi") %>%
  mutate(value = as.numeric(value)) %>%
  
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(pop1 ~ chromosome,
             scales = "free", switch = "x", space = "free_x") +
  xlab("Chromosome-scale scaffold") +
  ylab("pi") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.1, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        text = element_text(family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, by = 0.05),
                     expand = c(0, 0))

p1
ggsave("DretUKmerge_pixy_20260701_p1.png",
       units = "mm", width = 90, height = 90)

# -----------------------------
# WATTERSON THETA
# -----------------------------
p2 <- pixy_df %>%
  mutate(chrom_color_group = if_else(chromosome %% 2 == 1, "even", "odd")) %>%
  filter(statistic == "avg_watterson_theta") %>%
  mutate(value = as.numeric(value)) %>%
  
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(pop1 ~ chromosome,
             scales = "free", switch = "x", space = "free_x") +
  xlab("Chromosome-scale scaffold") +
  ylab("Watterson's theta") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.1, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        text = element_text(family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, by = 0.05),
                     expand = c(0, 0))

p2
ggsave("DretUKmerge_pixy_20260701_p2.png",
       units = "mm", width = 90, height = 90)

# -----------------------------
# TAJIMAS D
# -----------------------------
p3 <- pixy_df %>%
  mutate(chrom_color_group = if_else(chromosome %% 2 == 1, "even", "odd")) %>%
  filter(statistic == "tajima_d", pop1 == "Netherlands") %>%
  mutate(value = as.numeric(value)) %>%
  
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(pop1 ~ chromosome,
             scales = "free", switch = "x", space = "free_x") +
  xlab("Chromosome-scale scaffold") +
  ylab("Tajima's D") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.1, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        text = element_text(family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-7, 3),
                     breaks = seq(-7, 3, by = 2),
                     expand = c(0, 0))

p3
ggsave("DretUKmerge_pixy_20260701_p3.png",
       units = "mm", width = 90, height = 90)

# -----------------------------
# DXY
# -----------------------------
p4 <- pixy_df %>%
  mutate(chrom_color_group = if_else(chromosome %% 2 == 1, "even", "odd")) %>%
  filter(statistic == "avg_dxy") %>%
  mutate(value = as.numeric(value),
         pair = interaction(pop1, pop2, sep = " vs ")) %>%
  
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(pair ~ chromosome,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(value = label_value)) +
  xlab("Chromosome-scale scaffold") +
  ylab("Dxy") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.1, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        text = element_text(family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = seq(0, 0.25, by = 0.05),
                     # labels = c("0","0.3","0.6","0.9","0.8","1.1"),
                     expand = c(0, 0))

p4

ggsave("DretUKmerge_pixy_20260701_p4.png",
       units = "mm", width = 90, height = 90)

# -----------------------------
# FST HUDSON
# -----------------------------
p5 <- pixy_df %>%
  mutate(chrom_color_group = if_else(chromosome %% 2 == 1, "even", "odd")) %>%
  filter(statistic == "avg_hudson_fst") %>%
  mutate(value = as.numeric(value),
         pair = interaction(pop1, pop2, sep = " vs ")) %>%
  
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(pair ~ chromosome,
             scales = "free", switch = "x", space = "free_x",
             labeller = labeller(value = label_value)) +
  xlab("Chromosome-scale scaffold") +
  ylab("FST Hudson") +
  scale_color_manual(values = c("grey50", "black")) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.1, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        text = element_text(family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.6, 1.0),
                     breaks = c(-0.6, -0.4, -0.2, 0, 0.2,
                                0.4, 0.6, 0.8, 1.0),
                     expand = c(0, 0))

p5

ggsave("DretUKmerge_pixy_20260701_p5.png",
       units = "mm", width = 90, height = 90)

# -----------------------------
# ARRANGE ALL FIGS
# -----------------------------
figure <- ggarrange(p2,p3,p4,p5,
                    ncol = 2,
                    nrow = 2)

figure
ggsave("DretUKmerge_sup_v1.png", units = "mm", width = 180, height = 180)
ggsave("DretUKmerge_sup_v1.pdf", units = "mm", width = 180, height = 180, dpi = 600)

# -----------------------------
# OVERALL PI STATS
# -----------------------------
pixy_df %>%
  filter(statistic %in% c("count_diffs", "count_comparisons")) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1) %>%
  summarise(pi = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE)) %>%
  write_tsv("DretUKmerge_pi_summary.tsv")

pi_by_chrom <- pixy_df %>%
  filter(statistic %in% c("count_diffs", "count_comparisons")) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1, chromosome) %>%
  summarise(pi = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE),
            .groups = "drop")

write_tsv(pi_by_chrom, "DretUKmerge_pi_chromosome.tsv")

# -----------------------------
# OVERALL WATTERSON-THETA STATS
# -----------------------------
pixy_df %>%
  filter(statistic %in% c("raw_watterson_theta", "no_sites")) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value,
              values_fn = mean) %>%
  group_by(pop1) %>%
  summarise(wt = sum(raw_watterson_theta, na.rm = TRUE) / sum(no_sites, na.rm = TRUE)) %>%
  write_tsv("DretUKmerge_wt_summary.tsv")

wt_by_chrom <- pixy_df %>%
  filter(statistic %in% c("raw_watterson_theta", "no_sites")) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value,
              values_fn = mean) %>%
  group_by(pop1, chromosome) %>%
  summarise(wt = sum(raw_watterson_theta, na.rm = TRUE) / sum(no_sites, na.rm = TRUE),
            .groups = "drop")

write_tsv(wt_by_chrom, "DretUKmerge_wt_chromosome.tsv")

# -----------------------------
# OVERALL DXY STATS
# -----------------------------
pixy_df %>%
  filter(statistic %in% c("count_diffs", "count_comparisons"),
         !is.na(pop2)) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1, pop2) %>%
  summarise(dxy = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE),
            .groups = "drop") %>%
  write_tsv("DretUKmerge_dxy_summary.tsv")

dxy_by_chrom <- pixy_df %>%
  filter(statistic %in% c("count_diffs", "count_comparisons"),
         !is.na(pop2)) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1, pop2, chromosome) %>%
  summarise(dxy = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE),
            .groups = "drop")

write_tsv(dxy_by_chrom, "DretUKmerge_dxy_chromosome.tsv")

# -----------------------------
# OVERALL FST STATS
# -----------------------------
pixy_df %>%
  filter(statistic %in% c("hudson_fst_num", "hudson_fst_den"),
         !is.na(pop2)) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1, pop2) %>%
  summarise(fst = sum(hudson_fst_num, na.rm = TRUE) / sum(hudson_fst_den, na.rm = TRUE),
            .groups = "drop") %>%
  write_tsv("DretUKmerge_fst_summary.tsv")

fst_by_chrom <- pixy_df %>%
  filter(statistic %in% c("hudson_fst_num", "hudson_fst_den"),
         !is.na(pop2)) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  group_by(pop1, pop2, chromosome) %>%
  summarise(fst = sum(hudson_fst_num, na.rm = TRUE) / sum(hudson_fst_den, na.rm = TRUE),
            .groups = "drop")

write_tsv(fst_by_chrom, "DretUKmerge_fst_chromosome.tsv")

# -----------------------------
# OVERALL TAJIMA'S D STATS
# -----------------------------
# Load tajima_d files separately, keeping tajima_d_s_counts
tajima_files <- list.files(pixy_folder, pattern = "_tajima_d\\.txt$", full.names = TRUE)

tajima_raw <- map_dfr(tajima_files, ~ read_delim(.x, delim = "\t",
                                                 col_types = cols(tajima_d_s_counts = col_character(),
                                                                  .default = col_guess()))) %>%
  left_join(chromosome_key, by = c("chromosome" = "chr_name")) %>%
  mutate(chromosome = chr_numb) %>%
  select(-chr_numb)

parse_tajima_d_s_counts <- function(value) {
  if (length(value) == 0 || is.na(value)) {
    return(setNames(numeric(), character()))
  }
  value <- as.character(value)
  if (value == "" || value == "NA") {
    return(setNames(numeric(), character()))
  }
  counts <- setNames(numeric(), character())
  for (item in strsplit(value, ",", fixed = TRUE)[[1]]) {
    pair <- strsplit(item, ":", fixed = TRUE)[[1]]
    n <- pair[1]
    old <- counts[n]
    if (is.na(old)) old <- 0
    counts[n] <- old + as.numeric(pair[2])
  }
  counts
}

combine_tajima_d_s_counts <- function(values) {
  total <- setNames(numeric(), character())
  for (value in values) {
    counts <- parse_tajima_d_s_counts(value)
    for (n in names(counts)) {
      old <- total[n]
      if (is.na(old)) old <- 0
      total[n] <- old + counts[n]
    }
  }
  total
}

calc_tajima_d_stdev <- function(s_counts) {
  stdev <- 0
  for (n_name in names(s_counts)) {
    n <- as.integer(n_name)
    s <- as.numeric(s_counts[n_name])
    if (is.na(n) || n < 2 || s <= 0) next
    i <- seq_len(n - 1)
    a1 <- sum(1 / i)
    a2 <- sum(1 / (i^2))
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
    c1 <- b1 - (1 / a1)
    c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
    e1 <- c1 / a1
    e2 <- c2 / (a1^2 + a2)
    stdev <- stdev + sqrt((e1 * s) + (e2 * s * (s - 1)))
  }
  stdev
}

aggregate_tajima_d <- function(rows) {
  raw_pi <- sum(rows$raw_pi, na.rm = TRUE)
  raw_watterson_theta <- sum(rows$raw_watterson_theta, na.rm = TRUE)
  s_counts <- combine_tajima_d_s_counts(rows$tajima_d_s_counts)
  tajima_d_stdev <- calc_tajima_d_stdev(s_counts)
  tajima_d <- if (tajima_d_stdev <= 0) NA_real_ else (raw_pi - raw_watterson_theta) / tajima_d_stdev
  data.frame(tajima_d = tajima_d,
             no_sites = sum(rows$no_sites, na.rm = TRUE),
             raw_pi = raw_pi,
             raw_watterson_theta = raw_watterson_theta,
             tajima_d_stdev = tajima_d_stdev)
}

tajima_wide <- tajima_raw %>%
  filter(pop == "Netherlands") %>%
  select(pop, chromosome, window_pos_1, window_pos_2,
         raw_pi, raw_watterson_theta, no_sites, tajima_d_s_counts) %>%
  mutate(across(c(raw_pi, raw_watterson_theta, no_sites), as.numeric))

tajima_by_chrom <- tajima_wide %>%
  group_by(chromosome) %>%
  group_modify(~ aggregate_tajima_d(.x)) %>%
  ungroup()

write_tsv(tajima_by_chrom, "DretUKmerge_td_chromosome.tsv")

s_counts_genome <- combine_tajima_d_s_counts(tajima_wide$tajima_d_s_counts)
tajima_d_stdev_genome <- calc_tajima_d_stdev(s_counts_genome)
tajima_d_genome <- (sum(tajima_wide$raw_pi, na.rm = TRUE) - sum(tajima_wide$raw_watterson_theta, na.rm = TRUE)) / tajima_d_stdev_genome

data.frame(
  label = "genome_wide",
  tajima_d = tajima_d_genome,
  tajima_d_stdev = tajima_d_stdev_genome,
  raw_pi = sum(tajima_wide$raw_pi, na.rm = TRUE),
  raw_watterson_theta = sum(tajima_wide$raw_watterson_theta, na.rm = TRUE)
) %>%
  write_tsv("DretUKmerge_td_summary.tsv")
