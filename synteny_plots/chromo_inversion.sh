#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=tick_genespace_inversion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem=100G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skc49482@uga.edu

start_time=`date +%s`

module load GENESPACE/1.2.3-foss-2023a-R-4.3.2

## genespace
cat > genespace_inversions.R << 'EOF'
.libPaths(.libPaths()[!grepl("^/home/", .libPaths())])
library(GENESPACE)
library(ggplot2)
library(patchwork)

all_iso <- local({ load("/scratch/skc49482/Katie_genespace/genespace_aug/prep_input/genespace/results/gsParams.rda"); gsParam })
longest_iso <- local({ load("/scratch/skc49482/Katie_genespace/genespace_aug/prep_input/longest_isoform/genespace_isoform/results/gsParams.rda"); gsParam })

customPal <- colorRampPalette(c("#32A251FF", "#FF7F0FFF", "#3CB7CCFF",
                                "#FFD94AFF", "#39737CFF", "#B85A0DFF"))

build_plot <- function(gsParam, refG, gIDs, spLabels, flips, labFun = NULL, title = NULL) {
  args <- list(
    gsParam = gsParam, refGenome = refG, genomeIDs = gIDs,
    useRegions = FALSE,
    addThemes = theme(panel.background = element_rect(fill = "white")),
    braidAlpha = 0.75, palette = customPal, xlabel = "",
    syntenyWeight = 1, chrLabFontSize = 5,
    chrBorderCol = "black", chrBorderLwd = 0.1,
    invertTheseChrs = flips
  )
  if (!is.null(labFun)) args$chrLabFun <- labFun
  ripd <- do.call(plot_riparian, args)
  p <- ripd$plotData$ggplotObj
  glab <- unique(ripd$plotData$sourceData$chromosomes[, c("genome","y1","y2")])
  yb <- (glab$y1 + glab$y2)/2; names(yb) <- glab$genome
  pl <- parse(text = spLabels[names(yb)])
  pm <- p + scale_y_continuous(breaks = yb, labels = pl,
                               expand = c(0.01,0.01), name = NULL) +
    theme(axis.text.y = element_text(size = 6.5))
  if (length(pm$layers) >= 5) { pm$layers[[4]] <- NULL; pm$layers[[4]] <- NULL }
  if (!is.null(title)) pm <- pm + ggtitle(title) +
    theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))
  pm
}

flip_df <- data.frame(
  genome = c(
    "D_albipictus","D_albipictus","D_albipictus","D_albipictus",
    "D_andersoni","D_andersoni","D_andersoni","D_andersoni",
    "D_variabilis","D_variabilis","D_variabilis",
    "D_variabilis","D_variabilis","D_variabilis","D_variabilis",
    "D_silvarum","D_silvarum",
    "DretLouise","DretLouise",
    "DretUK","DretUK"
  ),
  chr = c(
    "9","10","5","8",
    "7","3","8","1",
    "7","10","11","5","9","3","6",
    "1","2",
    "1","2",
    "1","2"
  ),
  stringsAsFactors = FALSE
)

lab <- c(
  DretUK="italic('D. reticulatus')~'(Devon)'",
  DretLouise="italic('D. reticulatus')~'(Louise)'",
  D_silvarum="italic('D. silvarum')", D_albipictus="italic('D. albipictus')",
  D_andersoni="italic('D. andersoni')", D_variabilis="italic('D. variabilis')")
gids <- c("DretUK","DretLouise","D_silvarum","D_albipictus","D_andersoni","D_variabilis")

p_longest <- build_plot(longest_iso, "DretUK", gids, lab, flip_df)
p_alliso  <- build_plot(all_iso,     "DretUK", gids, lab, flip_df)

ggsave("Tick_longest_isoform_aug.pdf", p_longest, width=180, height=108, units="mm", dpi=600)
ggsave("Tick_alliso_aug.pdf",          p_alliso,  width=180, height=108, units="mm", dpi=600)
ggsave("Tick_combined_aug.pdf", p_alliso / p_longest, width=180, height=220, units="mm", dpi=600)

EOF
Rscript genespace_inversions.R

end_time=`date +%s`
runtime=$((end_time - start_time))
runtimeH=$((runtime / 3600))
runtimeM=$(((runtime % 3600) / 60))
runtimeS=$((runtime % 60))
echo "Duration: $runtime seconds"
echo "Duration: $runtimeH hours, $runtimeM minutes, $runtimeS seconds"
