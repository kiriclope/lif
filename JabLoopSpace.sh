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
	read dI CrecE CrecI CrecS CrecV Cff JabMin dJab JabMax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read dI CrecE CrecI CrecS Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read dI CrecE CrecI Cff JabMin dJab JabMax <<< "${6} ${7} ${8} ${9} ${10} ${11} ${12}"
    fi
    
    if [ "$nbpop" = 1 ]; then 
	read dI CrecI Cff dI Imax <<< "${6} ${7} ${8} ${9} ${10}"
    fi

fi

LANG=en_US

echo dI$dI JabMin$JabMin dJab$dJab JabMax$JabMax
for Jab in $(seq $JabMin $dJab $JabMax); do 
    
    DIRECTORY=$(printf "../Simulations/%dpop/%s/N%d/K%d/g%.2f/Prtr_I/Iext_I%.4f/Space/Crec4%.2fCrecI%.4fCff%.2f" "$nbpop" "$dir" "$N" "$k" "$g" "$I" "$CrecE" "$CrecI" "$Cff")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}CrecV${CrecV}Cff${Cff}Jab${Jab} ./ring2D.out $dir $nbpop $N $k $g $dI $CrecE $CrecI $CrecS $CrecV $Cff $Jab
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecE $CrecI $CrecS $Cff
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff}Jab${Jab} ./JabLoopSpace.out $dir $nbpop $N $k $g $dI $CrecE $CrecI $Cff $Jab 
	fi
	
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff} ./a.out $dir $nbpop $N $k $g $I $CrecI $Cff
	fi
	
	echo LIFspace_${dir}_N${N}K${k}Iext${I}CrecE${CrecE}CrecI${CrecI}Cff${Cff}Jab${Jab} "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
