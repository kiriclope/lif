#!/bin/bash
LANG=en_US
export LANG

if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N K g" dir nbpop N K g
else
    read dir nbpop N K g <<<  "$1 $2 $3 $4 $5"    
    read Imin dI Imax Prtr <<< "${6} ${7} ${8} ${9}"
fi


sed -i 's/IF_BENCHMARK 1/IF_BENCHMARK 0/' GlobalVars.h 
sed -i 's/IF_INTERPOLATION 1/IF_INTERPOLATION 0/' GlobalVars.h 
sed -i 's/IF_EULER 1/IF_EULER 0/' GlobalVars.h 
sed -i 's/IF_RK2 0/IF_RK2 1/' GlobalVars.h 

sed -i 's/IF_TRANSIENT_IEXT 1/IF_TRANSIENT_IEXT 0/' GlobalVars.h 
sed -i 's/IF_JabLoop 1/IF_JabLoop 0/' GlobalVars.h 
sed -i 's/IF_TIMECOURSE 1/IF_TIMECOURSE 0/' GlobalVars.h
sed -i 's/IF_SHARED 1/IF_SHARED 0/' GlobalVars.h 
 
sed -i 's/IF_IEXT 0/IF_IEXT 1/' GlobalVars.h 
sed -i 's/IF_GAUSS 1/IF_GAUSS 0/' GlobalVars.h 
sed -i 's/IF_RING 1/IF_RING 0/' GlobalVars.h 
sed -i 's/IF_SPEC 1/IF_SPEC 0/' GlobalVars.h 

sed -i 's/nbPref .*/nbPref 110971.5/' GlobalVars.h 
sed -i 's/IF_Nk 0/IF_Nk 1/' GlobalVars.h 
sed -i 's/DT .*/DT .01/' GlobalVars.h 
sed -i 's/IF_Prtr 0/IF_Prtr 1/' GlobalVars.h 
sed -i "s/PrtrPop .*/PrtrPop $Prtr/" GlobalVars.h 

sleep 1

echo "g++ ifnet.cpp -std=c++11 -Ofast -s -o IextLoop.out"
g++ ifnet.cpp -std=c++11 -Ofast -s -o IextLoop.out
sleep 5



echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Iext%s_%.4f/" "$nbpop" "$dir" "$N" "$K" "$g" "$Prtr" "$I") 
    
    if [ ! -d "$DIRECTORY" ]; then
        
	screen -dmS LIF_N${N}K${K}Iext${Prtr}_${I} ./IextLoop$.out $dir $nbpop $N $K $g $I 
	
	echo LIF_N${N}K${K}Iext${Prtr}_${I} $DIRECTORY "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
