#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh


#Input Files :- //////////////////////////////////////////////////////////////////////////////

id=$1
R1=$2 
R2=$3
cores=$4
output_report=$(echo "${id}_contamination_report.tsv")
bwadb=./bwadb/grh38_decoy/GRCh38_full_analysis_set_plus_decoy_hla.fa
krakendb=./krakendb/protocol_db



#Output Files :- /////////////////////////////////////////////////////////////////////////////

sam=./sam_file/BWA_${id}.sam
cram=./cram_file/BWA_${id}.cram
ubam=./ubam_file/BWA_${id}.bam
sorted_ubam=./ubam_file/sorted_BWA_${id}.qsort.bam
unampped_R1=./unmapped/unmapped_BWA_${id}_R1.fastq
unmapped_R2=./unmapped/unmapped_BWA_${id}_R2.fastq
kraken=./Kraken2/BWA_${id}.kraken2
kreport=./k2report/BWA_${id}.k2report
bracken=./bracken/BWA_${id}.bracken
breport=./breport/BWA_${id}.breport



#Help Desk :- ////////////////////////////////////////////////////////////////////////////////


if [ "$id" == "--help" ]; then   
    echo "script.sh id R1 R2 cores"
    exit 0
else




#Analysis 01 :- /////////////////////////////////////////////////////////////////////////////////
cd ./bwa; ./bwa mem -t $cores $bwadb $R1 $R2  > $sam; cd -
samtools view -@ $cores -T $bwadb  -C -o $cram $sam
samtools view -@ $cores -u -f12 -F256 $cram > $ubam
samtools sort -@ $cores -n -o $sorted_ubam $ubam
bedtools bamtofastq -i $sorted_ubam -fq $unampped_R1 -fq2 $unmapped_R2
kraken2 --db $krakendb --report $kreport --report-minimizer-data --minimum-hit-groups 3 $unmapped_R1 $unmapped_R2 > $kraken
bracken -d $krakendb -i $kreport -o $bracken -w $breport -r 100 -l S -t 10
rm -rf $sam

#Analysis 02 :- ///////////////////////////////////////////////////////////////////////////////


Total_read_count=$(less $R1 | awk 'END {print NR/4}')
Unmapped_count=$(less $unampped_R1 | awk ' END {print NR/4}')
Mapped_Read_count=$(awk -v t="$Total_read_count" -v u="$Unmapped_count" 'BEGIN {print t-u}')
Mapped_Read_coverage=$(awk -v m="$Mapped_Read_count" 'BEGIN {print (m*150)/3000000000}')
U_T=$(awk -v u="$Unmapped_count" -v t="$Total_read_count" 'BEGIN {print u/t}')
Total_coverage=$(awk -v t="$Total_read_count" 'BEGIN {print (t*150)/3000000000}')
microbe_count=$(less $breport | awk '$3 != 0   {print $6$7$8"\t"$3}' | grep -v Homosapiens | wc -l) 
Read_count=$(less $breport | awk ' $3 != 0  {print $6$7$8"\t"$3}' | grep -v "Homosapiens" | sort -k2,2n | tail -n 1 | awk '{print $2}')
Norm=$(awk -v m="$microbe_count" -v r="$Read_count" 'BEGIN {print r/m}')
Modest_microbe=$(less $breport | awk ' $3 != 0  {print $6$7$8"\t"$3}' | grep -v "Homosapiens" | sort -k2,2n | tail -n 1 | awk '{print $1}')

#Preprocessing :- ////////////////////////////////////////////////////////////////////////////////////////
if [ "$microbe_count" -eq 0 ]; then
    Norm=0
    Read_count=0
    Modest_microbe="NA"
fi

echo "Norm: $Norm"
echo "Read_count: $Read_count"
echo "Modest_microbe: $Modest_microbe"



#Analysis 03 :- //////////////////////////////////////////////////////////////////////////////////////////

conda run -n contaminon python3 ./python_contamination.py $id $Total_read_count $Unmapped_count $Mapped_Read_count $Mapped_Read_coverage $U_T $Total_coverage $Read_count $microbe_count $Norm $Modest_microbe $output_report 



fi
