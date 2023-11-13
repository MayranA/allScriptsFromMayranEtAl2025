
# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.7.1+galaxy0
# command_version:macs2 2.2.7.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'BAM filtered rmDup (as BED).bed'  --name 72h_WT_rep1    --format BED   --gsize '1870000000'      --SPMR     --call-summits  --keep-dup '1'  --d-min 20 --buffer-size 100000  --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp 72h_WT_rep1_peaks.xls 'individual macs2 tabular.tabular'   
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.7.1+galaxy0
# command_version:macs2 2.2.7.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'BAM filtered rmDup (as BED).bed'  --name 72h_WT_rep2    --format BED   --gsize '1870000000'      --SPMR     --call-summits  --keep-dup '1'  --d-min 20 --buffer-size 100000  --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp 72h_WT_rep2_peaks.xls 'individual macs2 tabular.tabular'   
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_intersectbed/2.30.0+galaxy1
# command_version:bedtools v2.30.0
bedtools intersect -a 'individual macs2 narrowPeaks.bed'  -b '72h_WT_rep2.bed'                 > 'Intersection of both peaks.bed'

# random_lines1
# command_version:
python '/data/galaxy/galaxy/server/tools/filters/random_lines_two_pass.py' 'BAM filtered rmDup (as BED).bed' 'Select random lines on data 119.bed' '49823996'

# random_lines1
# command_version:
python '/data/galaxy/galaxy/server/tools/filters/random_lines_two_pass.py' 'BAM filtered rmDup (as BED).bed' 'Select random lines on data 118.bed' '49823996'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.7.1+galaxy0
# command_version:macs2 2.2.7.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'Select random lines on data 118.bed' 'Select random lines on data 119.bed'  --name 72h_WT_rep1    --format BED   --gsize '1870000000'      --SPMR     --call-summits  --keep-dup '1'  --d-min 20 --buffer-size 100000   --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp 72h_WT_rep1_peaks.xls 'merged macs2 tabular.tabular'   
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_intersectbed/2.30.0+galaxy1
# command_version:bedtools v2.30.0
bedtools intersect -a 'merged macs2 narrowPeaks.bed'  -b 'Intersection of both peaks.bed'       -wa -wb          > 'bedtools Intersect intervals on data 229 and data 210.bed'

# Filter1
# command_version:
python '/data/galaxy/galaxy/server/tools/stats/filtering.py' 'bedtools Intersect intervals on data 229 and data 210.bed' 'Filter on data 230.bed' '/data/galaxy/galaxy/jobs/000/125/125889/configs/tmprx0444ze' 20 "str,int,int,str,int,str,float,float,float,int,str,int,int,str,int,str,float,float,float,int" 0

# Cut1
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/cutWrapper.pl' 'Filter on data 230.bed' 'c1-c10' T 'Cut on data 231.tabular'


# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sorted_uniq/1.1.0
# command_version:sort (GNU coreutils) 8.25
sort -u   -t '	' -o "narrowPeak shared by both replicates.bed" "Cut on data 231.tabular"
