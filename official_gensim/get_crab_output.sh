#!/bin/bash
echo "no merge list"
for D in crab*/
do
    
    cd ${D}
    echo "${D} beginning"

    scanWarning=1
    n=1
    while [ "${scanWarning}" == "1" ]
    do
        crab getoutput -d $PWD>../craboutput.txt
        grep "Warn" ../craboutput.txt || scanWarning=0
        echo "n = ${n}, scanWarning = $scanWarning"
        let n++
    done

    cd results
    #hadd NCUGlobalTuples.root *root
    echo "${D} done"
    #test -f NCUGlobalTuples.root && echo "exist" || echo $PATH
    # test the merge root file
    cd ../..
done
