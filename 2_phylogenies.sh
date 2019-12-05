#### make phylogenies of lamprologine ###
#pooja.singh09@gmail.com


#1 get fasta from bam files

sh bam2fasta.sh

#2 extract and align gene sequences

sh extractgeneseq.sh 

sh MSA.sh

#3 make raxml gene trees

sh raxml_genetrees.sh

#4 make astral species tree


