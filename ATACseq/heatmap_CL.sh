# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_intersectbed/2.30.0+galaxy1
# command_version:bedtools v2.30.0
bedtools intersect -a 'increased_96h.bed.bed'  -b 'narrowPeak shared by both replicates.bed'   -v   -wa          > 'increased_not_called_120h.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_intersectbed/2.30.0+galaxy1
# command_version:bedtools v2.30.0
bedtools intersect -a 'increased_96h.bed.bed'  -b 'narrowPeak shared by both replicates.bed'    -u  -wa          > 'increased_called_120h.bed'

# toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_compute_matrix/deeptools_compute_matrix/3.5.2+galaxy0
# command_version:computeMatrix 3.5.1
ln -f -s 'bigwigAverage on data 2391 and data 2390.bigwig' '48h_WT_0.bw' 
 ln -f -s 'bigwigAverage on data 2404.bigwig' '72h_WT_rep3_1.bw' 
 ln -f -s 'bigwigAverage on data 2407.bigwig' '96h_WT_rep3_2.bw' 
 ln -f -s 'bigwigAverage on data 2410.bigwig' '120h_WT_rep3_3.bw' 
 ln -f -s 'bigwigAverage on data 2406 and data 2405.bigwig' '72h_Snai1KO_4.bw' 
 ln -f -s 'bigwigAverage on data 2409 and data 2408.bigwig' '96h_Snai1KO_5.bw' 
 ln -f -s 'bigwigAverage on data 2411.bigwig' '120h_Snai1KO_Clone1_6.bw' 
    ln -f -s 'increased_called_120h.bed' 'increased_called_120h_0.bed' 
 ln -f -s 'increased_not_called_120h.bed' 'increased_not_called_120h_1.bed' 
  computeMatrix  reference-point --regionsFileName 'increased_called_120h_0.bed' 'increased_not_called_120h_1.bed'  --scoreFileName '48h_WT_0.bw' '72h_WT_rep3_1.bw' '96h_WT_rep3_2.bw' '120h_WT_rep3_3.bw' '72h_Snai1KO_4.bw' '96h_Snai1KO_5.bw' '120h_Snai1KO_Clone1_6.bw'  --outFileName 'computeMatrix on data 350, data 349, and others: Matrix.deeptools_compute_matrix_archive' --samplesLabel '48h_WT' '72h_WT_rep3' '96h_WT_rep3' '120h_WT_rep3' '72h_Snai1KO' '96h_Snai1KO' '120h_Snai1KO_Clone1'  --numberOfProcessors "${GALAXY_SLOTS:-4}"   --referencePoint center  --beforeRegionStartLength 2500 --afterRegionStartLength 2500  --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50     --transcriptID transcript --exonID exon --transcript_id_designator transcript_id

# toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_plot_heatmap/deeptools_plot_heatmap/3.5.2+galaxy0
# command_version:plotHeatmap 3.5.1
plotHeatmap --matrixFile 'computeMatrix on data 350, data 349, and others: Matrix.deeptools_compute_matrix_archive' --outFileName 'plotHeatmap on data 354: Image.eps'  --plotFileFormat 'eps'    --dpi '178'  --sortRegions 'no'   --sortUsing 'mean'  --averageTypeSummaryPlot 'mean'  --plotType 'lines'  --missingDataColor 'black'  --colorMap Greys  --alpha '1.0'   --sortUsingSamples 2  --xAxisLabel 'distance from peak center' --yAxisLabel 'genes'  --heatmapWidth 7.5 --heatmapHeight 25.0  --whatToShow 'heatmap and colorbar'  --startLabel 'TSS' --endLabel 'TES'  --refPointLabel 'peak center'     --legendLocation 'best'  --labelRotation '0'
