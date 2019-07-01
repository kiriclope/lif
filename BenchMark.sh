#!/usr/bin/env bash
LANG=en_US
export LANG

read dir nbpop N K g <<<  "$1 $2 $3 $4 $5" 

sed -i 's/IF_BENCHMARK 0/IF_BENCHMARK 1/' GlobalVars.h 
if [[ $dir == *"I"* ]] ; then
    sed -i 's/IF_INTERPOLATION 0/IF_INTERPOLATION 1/' GlobalVars.h 
else
    sed -i 's/IF_INTERPOLATION 1/IF_INTERPOLATION 0/' GlobalVars.h 
fi

if [[ $dir == *"EULER"* ]] ; then 
    echo $dir 
    sed -i 's/IF_EULER 0/IF_EULER 1/' GlobalVars.h 
    sed -i 's/IF_RK2 1/IF_RK2 0/' GlobalVars.h 
else
    echo $dir 
    sed -i 's/IF_EULER 1/IF_EULER 0/' GlobalVars.h 
    sed -i 's/IF_RK2 0/IF_RK2 1/' GlobalVars.h 	
fi

sed -i 's/IF_TRANSIENT_IEXT 1/IF_TRANSIENT_IEXT 0/' GlobalVars.h 
sed -i 's/IF_JabLoop 1/IF_JabLoop 0/' GlobalVars.h 
sed -i 's/IF_TIMECOURSE 1/IF_TIMECOURSE 0/' GlobalVars.h
sed -i 's/IF_SHARED 1/IF_SHARED 0/' GlobalVars.h 
 
sed -i 's/IF_GAUSS 1/IF_GAUSS 0/' GlobalVars.h 
sed -i 's/IF_RING 1/IF_RING 0/' GlobalVars.h 
sed -i 's/IF_SPEC 1/IF_SPEC 0/' GlobalVars.h 

sed -i 's/nbPref .*/nbPref 10000/' GlobalVars.h 
sed -i 's/IF_Nk 0/IF_Nk 1/' GlobalVars.h 
sed -i 's/DT .*/DT .01/' GlobalVars.h 
sed -i 's/IF_Prtr 1/IF_Prtr 0/' GlobalVars.h 

sleep 1

echo "g++ ifnet.cpp -std=c++11 -Ofast -s -o BenchMark.out"
g++ ifnet.cpp -std=c++11 -Ofast -s -o BenchMark.out
sleep 5

for dt in .001 .005 .01 .05 .1 .5 1; do 
    screen -dmS ${dir}_dt${dt} ./BenchMark.out $dir $nbpop $N $K $g $dt 
    echo ${dir}_dt${dt} "Running ..."
done

