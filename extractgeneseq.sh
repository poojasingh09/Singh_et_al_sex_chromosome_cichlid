for i in `ls *.fasta | cut -f 1 -d "."`;

do


/cl_tmp/singh/tools/samtools-1.3/bin/samtools faidx $i".fasta"

/cl_tmp/singh/tools/bedtools-2.25/bin/bedtools getfasta -name -fi $i".fasta" -bed Nbrichard.sexregion.bed -fo $i".candidate.fasta"

done
