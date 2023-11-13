path=/scratch/ldelisle/Alex_plots
pathRNA=${path}/toGEO/RNAseq
pathRNAsingle=${path}/toGEO/reanalysis_RNAseq/
path_FPKM_single_gastruloids=${path}/getData/single_gastruloids/AllCufflinks_Simplified.txt.gz
mkdir -p $path
cd $path
# Get gtf
wget "https://zenodo.org/record/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1" -O mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz -nc

# Generate bed for interesting genes:
zcat mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz | awk -F "\t" -v genes="Mixl1,Mesp1,Snai1,Cdh2,Cdh1,Eomes,T,Pou5f1,Pou3f1,Tbx6,Nodal" '
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
        width=end[gene]-start[gene]
        printf "%s\t%d\t%d\t%s\t0\t%s\n", chr[gene],start[gene] - width / 4, end[gene] + width / 4, gene, strand[gene]
        printf "%s\t%d\t%d\t%s\t0\t%s\n", chr[gene],start[gene], start[gene] + 1, gene, strand[gene] > "selected_genes_names.bed"
    }
}' > genes_to_plot.bed

# Single gastruloid sorted by Mesp1 expression:
zcat ${path_FPKM_single_gastruloids} | awk '
NR==1{
    for (i=4;i<=NF;i++) {
        gsub("FPKM_", "", $i)
        name[i] = $i
    }
}
$0~/Mesp1/{
    for (i=4;i<=NF;i++) {
        value[i] = $i
    }
}
END{
    for (i in name){
        print name[i]"\t"value[i]
    }
}' | sort -k2,2n | cut -f 1 > single_gastruloid_order.txt

ini_file=single_gastru.ini
levels=("" "d3" "c1" "b0" "9e" "8d" "7b" "6a" "58" "46" "35" "23")
echo "[scalebar]
where = top
fontsize = 20
file_type = scalebar
height = 1.5
" > ${ini_file}
for i in {1..10}; do 
    sample=$(awk -v i=$i NR==i'{print}' single_gastruloid_order.txt)
    level=${levels[$i]}
    echo "[$sample]
file = ${pathRNAsingle}/${sample}_both_strands.bw
color = #${level}${level}${level}
min_value = 0
max_value = TOFILL
height = 3
" >> ${ini_file}
done
echo "[genes]
file = mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
height = 1.5
color = black
border_color = black
prefered_name = gene_name
merge_transcripts = true
style = UCSC
merge_overlapping_exons = true
fontstyle = italic
arrow_interval = 20
arrowhead_fraction = 0.1
labels = false

[gene_names]
file = selected_genes_names.bed
fontstyle = italic
fontsize = 30
color = none
border_color = none
height = 3
" >> ${ini_file}

max_values=('' '40' '5')
for i in {1..2}; do
    gene=$(awk -v i=$i 'NR==i{print $4}' genes_to_plot.bed)
    region=$(awk -v i=$i 'NR==i{print $1":"$2"-"$3}' genes_to_plot.bed)
    sed "s/TOFILL/${max_values[$i]}/" ${ini_file} > ${ini_file/.ini/_${gene}.ini}
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.pdf} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.png} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.eps} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
done

## Mutants
for gene in Snai1 Cdh2; do
    ini_file=mutant_${gene}_positive.ini
    echo "[scalebar]
where = top
fontsize = 20
file_type = scalebar
height = 1.5
" > ${ini_file}
    if [ $gene = "Snai1" ]; then
        samples="96h_Snai1WT 96h_Snai1KO_Clone1 96h_Snai1KO_Clone2"
    elif [ $gene = "Cdh2" ]; then
        samples="120h_Cdh2Het 120h_Cdh2KO_Clone1 120h_Cdh2KO_Clone2"
    else
        echo "Unknown mutant"
        exit 1
    fi
    for sample in $samples; do
        if [[ "$sample" = *"WT"* ]]; then
            color=black
        elif [[ "$sample" = *"Het"* ]]; then
            color=black
        else
            color=red
        fi
        echo "[${sample}]
file = ${pathRNA}/${sample}_positive_strand.bw
title = ${sample} positive strand
color = ${color}
min_value = 0
max_value = TOFILL
height = 5
" >> ${ini_file}
    done
    echo "[genes]
file = mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
height = 0.9
color = black
border_color = black
prefered_name = gene_name
merge_transcripts = true
style = UCSC
display = collapsed
merge_overlapping_exons = true
arrow_interval = 20
arrowhead_fraction = 0.1
labels = false

[gene_names]
file = selected_genes_names.bed
fontstyle = italic
fontsize = 15
color = none
border_color = none
height = 1.5
" >> ${ini_file}
    cat ${ini_file} | sed 's/positive/negative/g' > ${ini_file/positive/negative}
done

genes=('Cdh2' 'Snai1')
max_values=('8' '10')
for i in 0 1; do
    gene=${genes[$i]}
    region=$(awk -v g=$gene '$4==g{print $1":"$2"-"$3}' genes_to_plot.bed)
    strand=$(awk -v g=$gene '$4==g{if($6=="+"){print "positive"}else{print "negative"}}' genes_to_plot.bed)
    ini_file=mutant_${gene}_${strand}.ini
    sed "s/TOFILL/${max_values[$i]}/" ${ini_file} > ${ini_file/.ini/_${gene}.ini}
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.pdf} --dpi 250 --plotWidth 5 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.png} --dpi 250 --plotWidth 5 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.eps} --dpi 250 --plotWidth 5 --fontSize 20 --trackLabelFraction 0.8
done

## Time-course
GreenLevels=('FF' 'E8' 'D2' 'BB' 'A5' '7B' '52' '29' '00')
ini_file=time_course_average_positive.ini
echo "[scalebar]
where = top
fontsize = 20
file_type = scalebar
height = 1.5
" > ${ini_file}
for i in {0..8}; do
    time=$((48+${i}*6))
    sample=${time}h_WT
    level=${GreenLevels[$i]}
    echo "[$sample]
file = ${pathRNA}/${sample}_positive_strand.bw
title = $sample
color = #FF${level}00
min_value = 0
max_value = TOFILL
height = 3
" >> ${ini_file}
done
echo "[genes]
file = mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
height = 1.5
color = black
border_color = black
prefered_name = gene_name
merge_transcripts = true
style = UCSC
merge_overlapping_exons = true
arrow_interval = 20
arrowhead_fraction = 0.1
labels = false
display = collapsed

[gene_names]
file = selected_genes_names.bed
fontstyle = italic
fontsize = 30
color = none
border_color = none
height = 3
" >> ${ini_file}
cat ${ini_file} | sed 's/positive/negative/g' > ${ini_file/positive/negative}

# Mixl1,Mesp1,Snai1,
# Cdh2,Cdh1,Eomes,
# T,Pou5f1,Pou3f1,
# Tbx6,Nodal
max_values=('' '' '' '12'
            '9' '40' '15'
            '200' '300' '30'
            '50' '30')
for i in {3..11}; do
    gene=$(awk -v i=$i 'NR==i{print $4}' genes_to_plot.bed)
    region=$(awk -v i=$i 'NR==i{print $1":"$2"-"$3}' genes_to_plot.bed)
    strand=$(awk -v g=$gene '$4==g{if($6=="+"){print "positive"}else{print "negative"}}' genes_to_plot.bed)
    ini_file=time_course_average_${strand}.ini
    sed "s/TOFILL/${max_values[$i]}/" ${ini_file} > ${ini_file/.ini/_${gene}.ini}
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.pdf} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.png} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
    pgt --tracks ${ini_file/.ini/_${gene}.ini} --region ${region} -o ${ini_file/.ini/_${gene}.eps} --dpi 250 --plotWidth 8 --fontSize 20 --trackLabelFraction 0.8
done
