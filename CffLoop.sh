#!/bin/bash


if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N k g" dir nbpop N k g
else
    read dir nbpop N k g <<<  "$1 $2 $3 $4 $5"

    if [ "$nbpop" = 4 ]; then 
	read I CrecE CrecI CrecS CrecV Cffmin dCff Cffmax Prtr <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read I CrecE CrecI CrecS Cffmin dCff Cffmax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read I CrecE CrecI Cffmin dCff Cffmax <<< "${6} ${7} ${8} ${9} ${10} ${11} "
    fi

    if [ "$nbpop" = 1 ]; then 
	read I CrecI Cffmin dCff Cffmax <<< "${6} ${7} ${8} ${9} ${10}"
    fi

fi

LANG=en_US

echo Cffmin$Cffmin dCff$dCff Cffmax$Cffmax
for Cff in $(seq $Cffmin $dCff $Cffmax); do
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Iext%s_%.4f/Space/Crec4%.2fCrecI%.4fCff%.2f" "$nbpop" "$dir" "$N" "$k" "$g" "$Prtr" "$I" "$CrecE" "$CrecI" "$Cff")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${Prtr}_${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}CrecV${CrecV}Cff${Cff} ./IextRing2D${Prtr}.out $dir $nbpop $N $k $g $I $CrecE $CrecI $CrecS $CrecV $Cff
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecE $CrecI $CrecS $Cff
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./d.out $dir $nbpop $N $k $g $I $CrecE $CrecI $Cff
	fi
	
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecI $Cff
	fi
	
	echo LIFspace_${dir}_N${N}K${k}Iext${Prtr}_${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
