path=/scratch/ldelisle/Alex_plots
pathATAC=/scratch/ldelisle/Alex_plots/toGEO/ATACseq
mkdir -p $path
cd $path
# Get gtf
wget "https://zenodo.org/record/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1" -O mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz -nc

# Get TSS for interesting genes:
zcat mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz | awk -F "\t" -v genes="Cdh1,Pou5f1,Nanog,Lefty1,Lefty2,Tcf15,Meox1,Foxc2" '
BEGIN{
    split(genes, a, ",")
    regex=""
    for (i in a) {
        start[a[i]] = 1000000000
        end[a[i]] = 0
        chr[a[i]] = "chr"
        if (regex != "") {
            regex = regex"|"a[i]
        } else {
            regex = a[i]
        }
    }
}
$9~"gene_name \""regex"\""{
    for (gene in start) {
        if ($9~"gene_name \""gene"\"") {
            chr[gene] = $1
            strand[gene] = $7
            if (start[gene] > $4) {
                start[gene] = $4
            }
            if (end[gene] < $5) {
                end[gene] = $5
            }
        }
    }
}
END{
    for (i in a) {
        gene = a[i]
        gene_strand = strand[gene]
        if (gene_strand == "+") {
            printf "%s\t%d\t%d\t%s\n", chr[gene],start[gene] - 1, start[gene], gene
        } else {
            printf "%s\t%d\t%d\t%s\n", chr[gene],end[gene] - 1, end[gene], gene
        }
    }
}' > genes_tss_to_plot.bed

# Manual for Esrrb
echo -e "chr12\t86420279\t86423599\tEsrrb" >> genes_tss_to_plot.bed

zcat mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz | grep "ENSMUST00000110203" > Esrrb.gtf

# Get peak overlapping TSS
bedtools intersect -a counts_on_peaks.txt.gz -b genes_tss_to_plot.bed -wa > ATAC_peak_tss.bed

# Increase 500bp each side to display
for f in ATAC_peak_tss.bed decreased_96h.bed increased_96h.bed; do
    bedtools slop -b 500 -i ${f} -g ~/genomes/fasta/mm10.fa.fai | bedtools sort | bedtools merge > ${f/.bed/_enlarged.bed}
done

ini_file=ATAC_Snai1.ini
echo "[scalebar]
where = top
fontsize = 20
file_type = scalebar
height = 1.5
" > ${ini_file}
for genotype in WT_rep3 Snai1KO; do
    for time in 72 96 120; do
        file=$(ls ${pathATAC}/${time}h_${genotype}*.bw | head -n 1)
        echo $file
        if [ ${genotype} = "Snai1WT" ]; then
            color="black"
        else
            color="red"
        fi
        echo "[$genotype $time]
file = $file
color = $color
title = $(basename $file .bw)
height = 3
summary_method = max
min_value = 0
max_value = TOFILL
" >> ${ini_file}
    done
done
echo "[genes]
file = mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
height = 0.5
color = black
border_color = black
prefered_name = gene_name
merge_transcripts = true
style = UCSC
merge_overlapping_exons = true
fontstyle = italic
arrow_interval = 20
arrowhead_fraction = 0.1
display = collapsed
labels = false

[gene_names]
file = genes_tss_to_plot.bed
fontstyle = italic
fontsize = 30
color = none
border_color = none
height = 3

[vhighlight1]
file = decreased_96h_enlarged.bed
type = vhighlight
color = blue
alpha = 0.3

[vhighlight2]
file = increased_96h_enlarged.bed
type = vhighlight
color = green
alpha = 0.3

[vhightlight3]
file = ATAC_peak_tss_enlarged.bed
type = vhighlight
color = red
alpha = 0.3
" >> ${ini_file}


genes=('Cdh1' 'Pou5f1' 'Nanog'
       'Tcf15' 'Meox1' 'Foxc2'
       'Lefty1-2' 'Esrrb')
regions=('chr8:106598000-106680000' 'chr17:35501000-35511000' 'chr6:122,699,521-122,760,831'
         'chr2:152115000-152155000' 'chr11:101888000-101898000' 'chr8:121110000-121385000'
         'chr1:180,864,109-180,965,689' 'chr12:86,391,493-86,556,012')
max_values=('20' '25' '20'
            '20' '32' '30'
            '25' '20')
plot_widths=('12' '4' '12'
             '6' '3' '18'
             '12' '12')

for i in {0..7}; do
    gene=${genes[$i]}
    region=${regions[$i]}
    plot_width=${plot_widths[$i]}
    sed "s/TOFILL/${max_values[$i]}/" ${ini_file} > ${ini_file/.ini/_${gene}.ini}
    if [ $gene = "Esrrb" ]; then
        sed -i "s/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz/Esrrb.gtf/" ${ini_file/.ini/_${gene}.ini}
    fi
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.pdf} --dpi 250 --plotWidth ${plot_width} --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.png} --dpi 250 --plotWidth ${plot_width} --fontSize 20 --trackLabelFraction 0.8
done
