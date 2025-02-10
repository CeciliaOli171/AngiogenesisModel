listrandomseed=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50) # list of random seed to test 
listsourceterm=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

for i in ${listrandomseed[@]}
do
    for j in ${listsourceterm[@]}
    do
        echo $i $j
        python3 /Users/coli171/Chaste/Source/projects/AngiogenesisModel/ModelAnalysis/ANZIAMPresentationTest3.py -valuerandomseed $i -sourceterm $j 
    done
done 