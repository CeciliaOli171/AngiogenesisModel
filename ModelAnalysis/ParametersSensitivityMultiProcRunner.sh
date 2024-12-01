listrandomseed=(11 22 33 44 55 66 77 88 99 100) # list of random seed to test 
listsigma=(2 3 4 5 6 7 8 9 10)
listomegap=(2 3 4 5 6 7 8 9 10)

for i in ${listrandomseed[@]}
do
    for j in ${listsigma[@]}
    do
        for k in ${listomegap[@]}
        do
            echo $i $j $k
            python3 /Users/coli171/Chaste/Source/projects/AngiogenesisModel/ModelAnalysis/ForcesCoeffTest.py -valuerandomseed $i -countersigma $j -counteromegap $k
        done
    done
done 