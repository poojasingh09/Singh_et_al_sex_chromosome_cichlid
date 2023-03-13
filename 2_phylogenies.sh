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


Singh et al 2023
_______________________________________________

ASTRAL: estimating species tree from gene tress
Mirabab et al 2014
________________________________________________

1. Concatenate all the IQtree gene trees
 cat *bestTree*

2. Run ASTRAL

java -Xmx10000M -jar /computation2/pooja/biotools/Astral/astral.4.7.12.jar -i all_raxml_bestTrees.tre -o all_IQtree_bestTrees_astralTree.tre

3. Run ASTRAL with 100 bootstraps


java -Xmx10000M -jar /computation2/pooja/biotools/Astral/astral.4.7.12.jar -i all_IQtree_bestTrees.tre -b path_bs -r 100 -o all_IQtree_bestTrees_astral.boot.site.out &

3.1 extract ASTRAL tree

tail -1 all_IQtree_bestTrees_astral.boot.site.out > all_IQtree_bestTrees_astral.boot.site.tre

4. Run ASTRAL with bootstraps (gene resampling): cant use 100 bootstraps: You
seem to have asked for 100 but only 71 replicaes could be created.
 Note that for gene resampling, you need more input bootstrap replicates than
the number of species tee replicates.


java -Xmx10000M -jar /computation2/pooja/biotools/Astral/astral.4.7.12.jar -i all_IQtree_bestTrees.tre -b path_bs -g -o all_IQtree_bestTrees_astral.boot.gene.out &

