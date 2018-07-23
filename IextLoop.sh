#!/bin/bash

S=0
V=0

if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N k g" dir nbpop N k g
    read -p "E I S V" E I S V
else
    read dir nbpop N k g <<<  "$1 $2 $3 $4 $5"
    
    if [ "$nbpop" = 4 ]; then 
	read E I S V dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read E I S dI Imax <<< "${6} ${7} ${8} ${9} ${10}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read E I dI Imax <<< "${6} ${7} ${8} ${9}"
    fi

    if [ "$nbpop" = 1 ]; then 
	read I dI Imax <<< "${6} ${7} ${8}"
    fi

fi

LANG=en_US

Imin=$I

echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Prtr_I/Iext_I%.4f/" "$nbpop" "$dir" "$N" "$k" "$g" "$I")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS LIF_${dir}_N${N}K${k}Iext${I} ./a.out $dir $nbpop $N $k $g $E $I $S $V
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS LIF_${dir}_N${N}K${k}Iext${I} ./a.out $dir $nbpop $N $k $g $E $I $S
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS LIF_N${N}K${k}Iext${I} ./a.out $dir $nbpop $N $k $g $E $I
	fi
	
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS LIF_${dir}_N${N}K${k}Iext${I} ./a.out $dir $nbpop $N $k $g $I
	fi
	
	echo LIF_N${N}K${k}Iext${I} $DIRECTORY "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
