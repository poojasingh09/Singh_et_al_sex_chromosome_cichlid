##pipeline for analysing lamprologini genome-resequencing data 
##pooja.singh09@gmail.com

for i in `cat a.txt | cut -d "." -f 1` ;  

do

echo "$i.starting" >> analysis_report.txt

######1. merge 8 bam files into one per sample

#qsub -pe mpi 8 -l h_vmem=5G -b y -cwd -N samMerge /usr/local/bioinf/bin/samtools merge /netapp/home/singh/pooja_phd/data/taborsky_callipterus/bourgeois_male_1.bam /nas2/projects/Neolamprologus/raw/24090*



######2. sort bam file

/usr/local/bioinf/bin/samtools sort -n /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.bam /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted

#rm $i.bam



#######3. extract fastq sequences from each file

/usr/local/bioinf/bin/bedtools bamtofastq -i /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted.bam -fq $i.sorted.1.fq -fq2 $i.sorted.2.fq

rm $i.sorted.bam

#bedtools bamtofastq -i bourgeois_male_1.bam -fq bourgeois_male_1.1.fq -fq2 bourgeois_male_1.2.fq




#######4 run fastQC quality control

##/usr/local/bioinf/bin/fastQC -o /netapp/home/singh/pooja_phd/data/taborsky_callipterus/ /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted.1.fq /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted.2.fq


######5 trim first 10 wonky reads

#qsub -pe mpi 8 -l h_vmem=5G -b y -cwd -N gtrim -o ./gtrim.log -e ./gtrim.err /bin/sh trim.sh 
#this is the script




#####5.1 trim first 10 bases which are wonky

/export/bio/fastx/latest/bin/fastx_trimmer -Q33 -f 9 -l 120 -i $i.sorted.1.fq -o $i.sorted.1.trimmed.fq

/export/bio/fastx/latest/bin/fastx_trimmer -Q33 -f 9 -l 120 -i $i.sorted.2.fq -o $i.sorted.2.trimmed.fq

rm $i.sorted.1.fq
rm $i.sorted.2.fq




#######5.2 get qual stats after trimming

##/export/bio/fastx/latest/bin/fastx_quality_stats -Q33 -i $i.sorted.1.trimmed.fq -o $i.sorted.1.trimmed.fq.stats

##/export/bio/fastx/latest/bin/fastx_quality_stats -Q33 -i $i.sorted.2.trimmed.fq -o $i.sorted.2.trimmed.fq.stats




######5.3 plot nucleotide distribution along read after trimming

##/export/bio/fastx/latest/bin/fastx_nucleotide_distribution_graph.sh -i $i.sorted.1.trimmed.fq.stats -o $i.sorted.1.trimmed.fq.statsND.png -t "fastq_nucleotide_distribution"

##/export/bio/fastx/latest/bin/fastx_nucleotide_distribution_graph.sh -i $i.sorted.2.trimmed.fq.stats -o $i.sorted.2.trimmed.fq.statsND.png -t "fastq_nucleotide_distribution"



########6 map to reference genome of N. brichardi 

#qsub -pe mpi 4 -l h_vmem=4G -b y -cwd -N cal_map -o ./cal_map.log -e ./cal_map.err /bin/sh map2ref.sh
#this is bwa+stampy; this method can be tried later; it is just so painfully SLOW

#qsub -pe mpi 8 -l h_vmem=5G -b y -cwd -V -N ON_map -o ./ON_map.log -e ./ON_map.err /bin/sh map2ref2.sh
#this is bowtie2 which is faster and apprently more accurate than bwa+stampy and thus will be used first 
##this is based on what was found with the heliconius genomes

/usr/local/bioinf/bin/bowtie2 --very-sensitive-local --no-mixed --no-discordant --threads 8 --time -x /netapp/home/singh/pooja_phd/data/cichlid_reference_genomes/Neolamprologous_brichardi/bt2_index/Nbri1.0 -1 $i.sorted.1.trimmed.fq -2 $i.sorted.2.trimmed.fq --un-gz $i.unpaired.unmapped.fq --un-conc-gz $i.discordant.unmapped.fq -S $i.Nb.bt2.sam

#rm $i.sorted.2.trimmed.fq
#rm $i.sorted.1.trimmed.fq



##########7 convert sam to bam 
/usr/local/bioinf/bin/samtools view -bS -T /netapp/home/singh/pooja_phd/data/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta $i.Nb.bt2.sam > $i.Nb.bt2.bam

#rm $i.Nb.bt2.sam




######8 sort bam
#qsub -pe mpi 8 -l h_vmem=5G -b y -cwd -N sortBAMfinal.$i -e sortBAMfinal.$i.err -o sortBAMfinal.$i.log 

#/usr/local/bioinf/bin/samtools sort $i.Nb.bt2.bam $i.Nb.bt2.sorted




#####9 flagstat on bam
#/usr/local/bioinf/bin/samtools flagstat $i.Nb.bt2.sorted.bam > $i.Nb.bt2.sorted.bam.flagstat


#10 add read groups

qsub -pe mpi 32 -b y -cwd -N RG -e rg.err -o rg.log /bin/sh read_group.sh 



#11 sort and index dedup bam

for i in `ls *RG.bam | cut -d "." -f -4` ;

do


        java -jar /cl_tmp/singh/tools/picard-tools-1.119/MarkDuplicates.jar I=$i.bam O=$i.dedup.bam M=$i.dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

	/cl_tmp/singh/tools/samtools-1.3/samtools sort $i.dedup.bam -o $i.dedup.sorted 

        #/cl_tmp/singh/tools/samtools-1.3/samtools index $i.dedup.sorted.bam

done


#12 indel realignment


for i in `ls *RG.dedup.sorted.bam | cut -d "." -f -6` ;

do


#echo $i.bam

#/cl_tmp/singh/tools/samtools-1.3/samtools index $i.bam


/usr/bin/java -jar -Xmx10G /cl_tmp/singh/tools/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta -I $i.bam -o $i.intervals -log $i.intervals.log

/usr/bin/java -jar -Xmx10G /cl_tmp/singh/tools/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta -I $i.bam -targetIntervals $i.intervals -o $i.realinged.bam -log $i.realign.log


#13 calcualte coverage

for i in `cat cal`;

do

/cl_tmp/singh/tools/bedtools-2.17.0/bin/bedtools coverage -abam $i -b /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/bedtools_taborsky_project/N_brichardi_v1.assembly.windows.30kb.s36.txt > $i".covD.s36.30kb.bed"

/cl_tmp/singh/tools/bedtools-2.17.0/bin/bedtools coverage -abam $i -b /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/bedtools_taborsky_project/N_brichardi_v1.assembly.windows.s23.30kb.15kbstep.txt > $i".covD.s23.30Kb.15kb.bed"

done



#14 GATK call SNPS and indels

for i in `ls *realinged.bam | cut -d "." -f -7` ;

do



/usr/bin/java -jar -Xmx100G /cl_tmp/singh/tools/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 16 -R /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta -I $i.bam --genotyping_mode DISCOVERY --minPruning 3 -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -stand_emit_conf 10 -stand_call_conf 30 -o $i.raw.vcf 

done


# filter SNPsc

for i in `ls *raw.vcf | cut -d '.' -f -7` ;

do 


#echo $i


/cl_tmp/singh/tools/vcftools_0.1.13/bin/vcftools --vcf $i.raw.vcf --out $i.filt --min-meanDP 8 --minQ 30 --remove-filtered-all --recode

done

## combine gVCF of filters vcfs

/usr/bin/java -Xms10G -Xmx320G -jar /cl_tmp/singh/tools/gatk/GenomeAnalysisTK.jar -T CombineGVCFs -R /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta --variant gVCF_list2.list -o combined.realigned.filt.recode.vcf

#15 make phlyogeny: this part is detailed in the next document


