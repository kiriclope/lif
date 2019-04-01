#!/bin/bash

if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N k g" dir nbpop N k g
    read -p "E I S V" E I S V
else
    read dir nbpop N k g <<<  "$1 $2 $3 $4 $5"

    if [ "$nbpop" = 4 ]; then 
	read Imin CrecE CrecI CrecS CrecV Cff dI Imax Prtr <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read Imin CrecE CrecI CrecS Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read Imin CrecE CrecI Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} "
    fi

    if [ "$nbpop" = 1 ]; then 
	read Imin CrecI Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10}"
    fi

fi

LANG=en_US

echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Iext%s_%.4f/Gauss2D/Crec4%.2fCrecI%.4fCff%.2f" "$nbpop" "$dir" "$N" "$k" "$g" "$Prtr" "$I" "$CrecE" "$CrecI" "$Cff")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${Prtr}_${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}CrecV${CrecV}Cff${Cff} ./IextRing2D${Prtr}.out $dir $nbpop $N $k $g $I $CrecE $CrecI $CrecS $CrecV $Cff
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecE $CrecI $CrecS $Cff
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./IextRing2D_EI.out $dir $nbpop $N $k $g $I $CrecE $CrecI $Cff
	fi
       
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecI $Cff
	fi
	
	echo LIFspace_${dir}_N${N}K${k}Iext${Prtr}_${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
