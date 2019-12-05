##pipeline for analysing lamprologini genome-resequencing data 
##pooja.singh09@gmail.com


for i in `cat bamfiles.txt | cut -d "." -f 1` ;  

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

/usr/local/bioinf/bin/fastQC -o /netapp/home/singh/pooja_phd/data/taborsky_callipterus/ /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted.1.fq /netapp/home/singh/pooja_phd/data/taborsky_callipterus/$i.sorted.2.fq


######5 trim first 10 wonky reads

qsub -pe mpi 8 -l h_vmem=5G -b y -cwd -N gtrim -o ./gtrim.log -e ./gtrim.err /bin/sh trim.sh 
#this is the script




#####5.1 trim first 10 bases which are wonky

/export/bio/fastx/latest/bin/fastx_trimmer -Q33 -f 9 -l 120 -i $i.sorted.1.fq -o $i.sorted.1.trimmed.fq

/export/bio/fastx/latest/bin/fastx_trimmer -Q33 -f 9 -l 120 -i $i.sorted.2.fq -o $i.sorted.2.trimmed.fq

rm $i.sorted.1.fq
rm $i.sorted.2.fq




#######5.2 get qual stats after trimming

/export/bio/fastx/latest/bin/fastx_quality_stats -Q33 -i $i.sorted.1.trimmed.fq -o $i.sorted.1.trimmed.fq.stats

/export/bio/fastx/latest/bin/fastx_quality_stats -Q33 -i $i.sorted.2.trimmed.fq -o $i.sorted.2.trimmed.fq.stats




######5.3 plot nucleotide distribution along read after trimming

/export/bio/fastx/latest/bin/fastx_nucleotide_distribution_graph.sh -i $i.sorted.1.trimmed.fq.stats -o $i.sorted.1.trimmed.fq.statsND.png -t "fastq_nucleotide_distribution"

/export/bio/fastx/latest/bin/fastx_nucleotide_distribution_graph.sh -i $i.sorted.2.trimmed.fq.stats -o $i.sorted.2.trimmed.fq.statsND.png -t "fastq_nucleotide_distribution"



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

rm $i.Nb.bt2.sam




######8 sort bam

/usr/local/bioinf/bin/samtools sort $i.Nb.bt2.bam $i.Nb.bt2.sorted




#####9 flagstat on bam
#/usr/local/bioinf/bin/samtools flagstat $i.Nb.bt2.sorted.bam > $i.Nb.bt2.sorted.bam.flagstat


#10 add read groups

qsub -pe mpi 32 -b y -cwd -N RG -e rg.err -o rg.log /bin/sh read_group.sh 



#11 sort and index dedup bam


sh dedup.sh



#12 samtools call SNPS and indels


sh varcall_all_Onil2.sh

echo "$i.finished" >> analysis_report.txt

done
