#!/bin/bash
echo "no merge list"
for D in */
do
    # echo "$D"
    cd ${D}
    crab getoutput -d $PWD
    cd results
    hadd NCUGlobalTuples.root *root
    #test -f NCUGlobalTuples.root && echo "exist" || echo $PATH
    # test the merge root file
    cd ../..
done
