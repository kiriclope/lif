#!/bin/bash

if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N K g" dir nbpop N K g
else
    read dir nbpop N K g <<<  "$1 $2 $3 $4 $5"    
    read Imin dI Imax Prtr <<< "${6} ${7} ${8} ${9}"
fi

LANG=en_US

echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Iext%s_%.4f/" "$nbpop" "$dir" "$N" "$k" "$g" "$Prtr" "$I") 
    
    if [ ! -d "$DIRECTORY" ]; then
        
	screen -dmS LIF_${dir}_N${N}K${K}Iext${Prtr}_${I} ./IextLoop${Prtr}.out $dir $nbpop $N $K $g $I 
	
	echo LIF_N${N}K${K}Iext${Prtr}_${I} $DIRECTORY "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
