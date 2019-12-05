for i in `ls *.fas`;

do

#echo $i
/cl_tmp/singh/tools/raxml/raxmlHPC-PTHREADS -T 10 -f a -s $i -p 589375 -x 12345 -N 100 -m GTRCAT -w /cl_tmp/singh/taborsky_callipterus/mapped_data/raxml_phylogeny/aligned/gene_trees/ -n $i

done
