#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=tick_genespace_longestIsoform
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

wd="/scratch/skc49482/Katie_genespace/genespace_aug/prep_input/longest_isoform/genespace_isoform"  

# running orthofinder
orthofinder -f $wd/peptide -t 16 -a 1 -X -o $wd/orthofinder

# run genespace
cat > run_gs.R << 'EOF'
.libPaths(.libPaths()[!grepl("^/home/", .libPaths())])
library(GENESPACE); library(ggplot2)
wd <- "/scratch/skc49482/Katie_genespace/genespace_aug/prep_input/longest_isoform/genespace_isoform"
path2mcscanx <- "/apps/eb/MCScanX/1.0.0-GCC-12.3.0/"

gpar <- init_genespace(wd=wd, path2mcscanx=path2mcscanx)
out  <- run_genespace(gpar, overwrite=TRUE)

customPal <- colorRampPalette(c("#32A251FF","#FF7F0FFF","#3CB7CCFF","#FFD94AFF","#39737CFF","#B85A0DFF"))


ripd <- plot_riparian(
  gsParam=out, refGenome="DretUK",
  genomeIDs=c("DretUK","DretLouise","D_silvarum","D_albipictus","D_andersoni","D_variabilis"),
  useRegions=FALSE, syntenyWeight=1, braidAlpha=0.75, palette=customPal, xlabel="",
  chrLabFontSize=5, chrBorderCol="black", chrBorderLwd=0.1,
  addThemes=theme(panel.background=element_rect(fill="white")))

p <- ripd$plotData$ggplotObj
glab <- unique(ripd$plotData$sourceData$chromosomes[,c("genome","y1","y2")])
yb <- (glab$y1+glab$y2)/2; names(yb) <- glab$genome
sp_labels <- c(DretUK="italic('D. reticulatus')~'(Devon)'",
  DretLouise="italic('D. reticulatus')~'(Louise)'",
  D_silvarum="italic('D. silvarum')", D_albipictus="italic('D. albipictus')",
  D_andersoni="italic('D. andersoni')", D_variabilis="italic('D. variabilis')")
p <- p + scale_y_continuous(breaks=yb, labels=parse(text=sp_labels[names(yb)]),
        expand=c(0.01,0.01), name=NULL) + theme(axis.text.y=element_text(size=6.5))

ggsave(file.path(wd,"riparian_longestiso.pdf"), p, width=180, height=108, units="mm", dpi=600)
EOF
Rscript run_gs.R

end_time=`date +%s`
runtime=$((end_time - start_time))
runtimeH=$((runtime / 3600))
runtimeM=$(((runtime % 3600) / 60))
runtimeS=$((runtime % 60))
echo "Duration: $runtime seconds"
echo "Duration: $runtimeH hours, $runtimeM minutes, $runtimeS seconds"
