# cat1
# command_version:
python /data/galaxy/galaxy/server/tools/filters/catWrapper.py 'Concatenate datasets on data 233, data 240, and others.bed' 'narrowPeak shared by both replicates.bed' 'narrowPeak shared by both replicates.bed' 'narrowPeak shared by both replicates.bed' 'narrowPeak shared by both replicates.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
sortBed -i 'Concatenate datasets on data 233, data 240, and others.bed'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'SortBed on Concatenate datasets on data 233, data 240, and others.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_mergebed/2.30.0
# command_version:bedtools v2.30.0
mergeBed -i 'SortBed on Concatenate datasets on data 233, data 240, and others.bed'  -d 0    > 'Merged SortBed on Concatenate datasets on data 233, data 240, and others.bed'

# toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_multi_bam_summary/deeptools_multi_bam_summary/3.5.2+galaxy0
# command_version:multiBamSummary 3.5.1
ln -s 'BAM filtered rmDup.bam' './0.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17235.dat' './0.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './1.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17236.dat' './1.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './2.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17237.dat' './2.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './3.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17238.dat' './3.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './4.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17239.dat' './4.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './5.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17240.dat' './5.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './6.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17241.dat' './6.bam.bai' 
 ln -s 'BAM filtered rmDup.bam' './7.bam' 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17242.dat' './7.bam.bai' 
   multiBamSummary BED-file --numberOfProcessors "${GALAXY_SLOTS:-4}"  --outFileName 'multiBamSummary on data 260, data 8, and others: correlation matrix.deeptools_coverage_matrix' --bamfiles '0.bam' '1.bam' '2.bam' '3.bam' '4.bam' '5.bam' '6.bam' '7.bam' --labels '129SV_48h_rep1' '129SV_48h_rep2' '129SV_72h_rep1' '129SV_72h_rep2' '129SV_96h_rep1' '129SV_96h_rep2' '129SV_120h_rep1' '129SV_120h_rep2'  --outRawCounts 'multiBamSummary on data 260, data 8, and others: bin counts.tabular'   --BED Merged SortBed on Concatenate datasets on data 233, data 240, and others.bed
