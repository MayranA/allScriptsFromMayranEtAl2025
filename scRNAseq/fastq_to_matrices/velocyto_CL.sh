
# toolshed.g2.bx.psu.edu/repos/iuc/velocyto_cli/velocyto_cli/0.17.17+galaxy1
# command_version: velocyto, version 0.17.17
mkdir -p '72h_WT_rep1/outs/filtered_gene_bc_matrices/whatever/' 
 ln -s 'RNA STARSolo on data 47, data 825, and others: Alignments.bam' '72h_WT_rep1/outs/possorted_genome_bam.bam' 
 ln -s 'DropletUtils 10X Barcodes on data 49, data 48, and data 50.tsv' '72h_WT_rep1/outs/filtered_gene_bc_matrices/whatever/barcodes.tsv' 
 velocyto  run10x  -t 'uint16'  --samtools-threads ${GALAXY_SLOTS:-1} --samtools-memory ${GALAXY_MEMORY_MB:-100}  '-vv' '72h_WT_rep1' 'bf5a119_mm10_allGastruloids_min10_extended.gtf.gz.gtf' 
 mv '72h_WT_rep1/velocyto/'*.loom 'output.loom'
