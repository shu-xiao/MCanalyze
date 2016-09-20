#!/bin/bash

filename='root_file.txt'
ls -d */>>$filename
exec< $filename
#line='crab_MonoHbb_ZpBaryonic_MZp-300_MChi-1_13TeV/results/NCUGlobalTuples.root'
#echo $line
while read line
do
    echo "${line}results/NCUGlobalTuples.root"
    root -q -b -l -n resolved_xAna_monoHiggs.C++\(\"${line}results/NCUGlobalTuples.root\"\)>test1.txt
    root -q -b -l -n boosted_xAna_monoHiggs.C++\(\"${line}results/NCUGlobalTuples.root\"\)>test2.txt
    echo 'transfering ...'
    sed '1,2d' test1.txt>"${line%/}_resolved.txt"
    sed '1,2d' test2.txt>"${line%/}_boosted.txt"
    echo "${line} done"    
done
rm test*.txt
echo 'finished'
