for i in `ls *.fa | cut -f -4 -d "."`;

do

/cl_tmp/singh/tools/mafft-7.273-with-extensions/bin/mafft $i".fa" > $i".aligned.fa"

/usr/bin/java -jar /cl_tmp/singh/tools/BMGE-1.12/BMGE.jar -i $i".aligned.fa" -t DNA -of $i".aligned.faa"

/usr/bin/java -jar /cl_tmp/singh/tools/BMGE-1.12/BMGE.jar -i $i".aligned.fa" -t DNA -op $i".aligned.phylip"

done
