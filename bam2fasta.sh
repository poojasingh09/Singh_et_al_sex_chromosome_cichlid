for i in `ls *.Nb.bt2.RG.dedup.sorted.realinged.bam`; 

do

echo qsub -l h_vmem=100G -b y -cwd -N bam2fasta -e b2f.err -o b2f.log /cl_tmp/singh/tools/angsd/angsd -doFasta 2 -doCounts 1 -i $i -out $i".fasta"

done
