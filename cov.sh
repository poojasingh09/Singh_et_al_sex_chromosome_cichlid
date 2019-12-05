for i in `ls *Nb.bt2.RG.dedup.sorted.realinged.bam`;

do


qsub -q all.q -l h_vmem=20G -b y -cwd -N cov -e cov.err -o cov.log /cl_tmp/singh/tools/bedtools-2.17.0/bin/bedtools coverage -abam $i -b /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/bedtools_taborsky_project/N_brichardi_v1.assembly.windows.s23.100kb.50kbstep.txt > $i".covD.s23.100Kb.50kb.bed"

/cl_tmp/singh/tools/bedtools-2.17.0/bin/bedtools coverage -abam ../../dwarf_male_1_daughter.Nb.bt2.RG.dedup.sorted.realinged.bam -b /cl_tmp/singh/cichlid_reference_genomes/Neolamprologous_brichardi/bedtools_taborsky_project/N_brichardi_v1.assembly.windows.1kb.txt > ../../dwarf_male_1_daughter.Nb.bt2.RG.dedup.sorted.realinged.bam.covD.1Kb.s23.bed


done

