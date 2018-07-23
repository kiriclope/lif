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
	read E I S V CrecE CrecI CrecS CrecV Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15} ${16}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read E I S CrecE CrecI CrecS Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read E I CrecE CrecI Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12}"
    fi

    if [ "$nbpop" = 1 ]; then 
	read I CrecI Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10}"
    fi

fi

LANG=en_US

Imin=$I

echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Prtr_I/Iext_I%.4f/Space/Crec4%.2fCrecI%.4fCff%.2f" "$nbpop" "$dir" "$N" "$k" "$g" "$I" "$CrecE" "$CrecI" "$Cff")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}CrecV${CrecV}Cff${Cff} ./a.out $dir $nbpop $N $k $g $E $I $S $V $CrecE $CrecI $CrecS $CrecV $Cff
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}Cff${Cff} ./a.out $dir $nbpop $N $k $g $E $I $S $CrecE $CrecI $CrecS $Cff
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./a.out $dir $nbpop $N $k $g $E $I $CrecE $CrecI $Cff
	fi
	
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecI $Cff
	fi
	
	echo LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
