#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=tick_genespace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem=400G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skc49482@uga.edu

start_time=`date +%s`

module load GENESPACE/1.2.3-foss-2023a-R-4.3.2

## input files
input_dir="/scratch/skc49482/Katie_genespace/Tick"
peptide="$input_dir/peptide"
orthofinder="$input_dir/orthofinder"

## OrthoFinder

echo -e "\n Running orthofinder \n"
orthofinder -f $peptide -t 16 -a 1 -X -o $orthofinder

echo -e "\n running genespace \n"

## genespace
cat > genespace_run_Tick.R << 'EOF'
.libPaths(.libPaths()[!grepl("^/home/", .libPaths())])
library(GENESPACE)
library(ggplot2)

## define wd here again
wd <- "/scratch/skc49482/Katie_genespace/Tick"
path2mcscanx <- "/apps/eb/MCScanX/1.0.0-GCC-12.3.0/"
gpar <- init_genespace(wd = wd, path2mcscanx = path2mcscanx)
out <- run_genespace(gpar, overwrite = T)

customPal <- colorRampPalette(c("#32A251FF", "#FF7F0FFF", "#3CB7CCFF",
                                "#FFD94AFF", "#39737CFF", "#B85A0DFF"))
                                
ggthemes <- theme(panel.background = element_rect(fill = "white"))
ripd <- plot_riparian(gsParam=out,
        refGenome="DretUK", 
        genomeIDs= c("DretUK","DretElska","Dsil",
                    "Dalb","Dand","Dvar"),
        useRegions= FALSE,
        addThemes  = ggthemes,
        braidAlpha = 0.75,
        palette    = customPal,
         xlabel     = "",
         invertTheseChrs = data.frame(
    genome = c(
      "Dalb", "Dand",
      "Dvar", "Dvar", "Dvar", "Dvar",
      "Dvar", "Dalb", "Dalb", "Dand",
      "Dvar", "Dvar", "Dalb", "Dalb",
      "Dalb", "Dand", "Dvar",
      "Dand", "Dvar"
    ),
    chr = c(
      "1", "7",
      "7", "10", "11", "5",
      "9", "9", "10", "3",
      "1", "3", "9", "5",
      "8", "8", "6",
      "2", "2"
    ),
    stringsAsFactors = FALSE
  )
)
# Italic labels + custom scale bar
p <- ripd$plotData$ggplotObj
glab <- unique(ripd$plotData$sourceData$chromosomes[, c("genome", "y1", "y2")])
y_breaks <- (glab$y1 + glab$y2) / 2
names(y_breaks) <- glab$genome
sp_labels <- c(
  "DretUK"       = "italic('D. reticulatus')~'(UK)'",
  "DretElska"       = "italic('D. reticulatus')~'(Elska)'",
  "Dsil"    = "italic('D. silvarum')",
  "Dalb"  = "italic('D. albipictus')",
  "Dand"   = "italic('D. andersoni')",
  "Dvar"  = "italic('D. variabilis')"
)
parsed_labels <- parse(text = sp_labels[names(y_breaks)])
p_mod <- p +
  scale_y_continuous(
    breaks = y_breaks,
    labels = parsed_labels,
    expand = c(0.01, 0.01),
    name   = NULL
  ) +
  theme(axis.text.y = element_text(size = 11))
p_mod2 <- p_mod
p_mod2$layers[[4]] <- NULL
p_mod2$layers[[4]] <- NULL

ggsave(
  filename = "Tick_riparian_plot.pdf",
  plot     = p_mod2,
  device   = "pdf",
  width    = 180,
  height   = 108,
  units    = "mm",
  dpi      = 600
)
EOF
Rscript genespace_run_Tick.R

end_time=`date +%s`
runtime=$((end_time - start_time))
runtimeH=$((runtime / 3600))
runtimeM=$(((runtime % 3600) / 60))
runtimeS=$((runtime % 60))
echo "Duration: $runtime seconds"
echo "Duration: $runtimeH hours, $runtimeM minutes, $runtimeS seconds"
