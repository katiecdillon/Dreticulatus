#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=tick_genespace_isoform
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
input_dir="/scratch/skc49482/Katie_genespace/Tick_isoform/genespace"
peptide="$input_dir/peptide"
orthofinder="$input_dir/orthofinder"

## orthofinder
echo -e "\n Running orthofinder \n"
orthofinder -f $peptide -t 16 -a 1 -X -o $orthofinder

echo -e "\n running genespace \n"

## genespace
cat > genespace_run_Tick_isoform.R << 'EOF'
.libPaths(.libPaths()[!grepl("^/home/", .libPaths())])
library(GENESPACE)
library(ggplot2)

## define wd here again
wd <- "/scratch/skc49482/Katie_genespace/Tick_isoform/genespace"
path2mcscanx <- "/apps/eb/MCScanX/1.0.0-GCC-12.3.0/"

original_fn <- get("run_mcscanx", envir = getNamespace("GENESPACE"))
fn_text <- deparse(original_fn)
fn_text <- gsub(
  'suppressWarnings(collin <- fread(cmd = sprintf("cat %s | grep %s_ | grep :",',
  'suppressWarnings(collin <- tryCatch(fread(cmd = sprintf("cat %s | grep %s_ | grep :",',
  fn_text, fixed = TRUE
)
fn_text <- gsub(
  'select = 1:3, showProgress = FALSE, header = FALSE))',
  'select = 1:3, showProgress = FALSE, header = FALSE), error = function(e) data.table::data.table(blkID = character(0), gn1 = character(0), gn2 = character(0))))',
  fn_text, fixed = TRUE
)
patched_fn <- eval(parse(text = fn_text))
environment(patched_fn) <- getNamespace("GENESPACE")
assignInNamespace("run_mcscanx", patched_fn, ns = "GENESPACE")

gpar <- init_genespace(wd = wd, path2mcscanx = path2mcscanx)
out <- run_genespace(gpar, overwrite = T)

customPal <- colorRampPalette(c("#32A251FF", "#FF7F0FFF", "#3CB7CCFF",
                                "#FFD94AFF", "#39737CFF", "#B85A0DFF"))
                                
ggthemes <- theme(panel.background = element_rect(fill = "white"))

ripd <- plot_riparian(gsParam=out,
        refGenome="Dret_UK", 
        genomeIDs= c("Dret_UK","Dret_EU","D_silvarum",
                                 "D_albipictus","D_andersoni","D_variabilis"),
        useRegions= FALSE,
        addThemes  = ggthemes,
        braidAlpha = 0.75,
        palette    = customPal,
         xlabel     = "",
         invertTheseChrs = data.frame(
    genome = c(
      "D_albipictus", "D_andersoni",
      "D_variabilis", "D_variabilis", "D_variabilis", "D_variabilis",
      "D_variabilis", "D_albipictus", "D_albipictus", "D_andersoni",
      "D_variabilis", "D_variabilis", "D_albipictus", "D_albipictus",
      "D_albipictus", "D_andersoni", "D_variabilis",
      "D_andersoni", "D_variabilis"
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
  "Dret_UK"       = "italic('D. reticulatus')~'(UK)'",
  "Dret_EU"       = "italic('D. reticulatus')~'(Elska)'",
  "D_silvarum"    = "italic('D. silvarum')",
  "D_albipictus"  = "italic('D. albipictus')",
  "D_andersoni"   = "italic('D. andersoni')",
  "D_variabilis"  = "italic('D. variabilis')"
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
  filename = "Tick_riparian_isoform_plot.pdf",
  plot     = p_mod2,
  device   = "pdf",
  width    = 180,
  height   = 108,
  units    = "mm",
  dpi      = 600
)

EOF
Rscript genespace_run_Tick_isoform.R

end_time=`date +%s`
runtime=$((end_time - start_time))
runtimeH=$((runtime / 3600))
runtimeM=$(((runtime % 3600) / 60))
runtimeS=$((runtime % 60))

echo "Duration: $runtime seconds"
echo "Duration: $runtimeH hours, $runtimeM minutes, $runtimeS seconds"
