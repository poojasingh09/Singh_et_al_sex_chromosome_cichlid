for i in `ls *.fa`;

do

sed -i -e 's/\(>\).*\(Nbri\)/\1\2/' $i 
sed -i -e 's/\(>\).*\(bour\)/\1\2/' $i 
sed -i -e 's/\(>\).*\(alto\)/\1\2/' $i
sed -i -e 's/\(>\).*\(lema\)/\1\2/' $i
sed -i -e 's/\(>\).*\(multif\)/\1\2/' $i 
sed -i -e 's/\(>\).*\(ocellatus\)/\1\2/' $i
sed -i -e 's/\(>\).*\(brevis\)/\1\2/' $i

done
