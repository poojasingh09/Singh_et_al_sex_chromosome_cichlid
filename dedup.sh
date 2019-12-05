for i in `ls *RG.bam | cut -d "." -f -4` ;

do


        java -jar /cl_tmp/singh/tools/picard-tools-1.119/MarkDuplicates.jar I=$i.bam O=$i.dedup.bam M=$i.dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

	/cl_tmp/singh/tools/samtools-1.3/samtools sort $i.dedup.bam -o $i.dedup.sorted 

        /cl_tmp/singh/tools/samtools-1.3/samtools index $i.dedup.sorted.bam

done
