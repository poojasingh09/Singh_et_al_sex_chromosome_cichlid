for i in `cat a`;

do

	while read line
	do
    	if [[ ${line:0:1} == '>' ]]
    	then
        	outfile=${line#>}.fa
        	echo $line"_$i" >> $outfile
    	else
        	echo $line >> $outfile
    	fi
	done < $i

done
