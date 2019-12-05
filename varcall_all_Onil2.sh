/cl_tmp/singh/tools/samtools-1.3/bin/samtools mpileup -q 30 -Q 30 -E -C 50 -t DP -t SP -m 40 -F 0.003164557 -u -f /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/N_brichardi_v1.assembly.fasta -b allbam.txt | /cl_tmp/singh/tools/bcftools-1.3/bin/bcftools call -p 0.999089 -vm - > all.lamprologini.raw2.bcf

/cl_tmp/singh/tools/bcftools-1.3/bin/bcftools view all.callip.raw.bcf | /cl_tmp/singh/tools/bcftools-1.3/bin/vcfutils.pl varFilter -d 8 -D 100000 -w0 -W0 -10 -20 -40 -e0 > ./all.callip.flt.vcf
