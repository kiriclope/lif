#!/bin/bash

S=0
V=0

if [ -z "$1" ]; then
    echo "Network"
    read -p "dir nbpop N" dir nbpop N
    read -p "dt duration Tst Tw g" dt duration Tst Tw g
    read -p "E I S V" E I S V
else
    read dir nbpop N k dt duration Tst Tw g <<<  "$1 $2 $3 $4 $5 $6 $7 $8 $9"

    if [ "$nbpop" = 4 ]; then 
	read E I S V  Tme Tmi Tms Tmv dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read E I S Tme Tmi Tms dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read E I Tme Tmi dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15}"
    fi

fi

export $LANG=en_US

Imin=$I 

echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do

    DIRECTORY=$(printf "../%dpop/Network/%s/N%d/K%d/g%.2f/Prtr_I/Iext_I%.4f/" "$nbpop" "$dir" "$N" "$k" "$g" "$I")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./ifnetIext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $S $V $Tme $Tmi $Tms $Tmv
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./ifnet_Iext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $S $Tme $Tmi $Tms
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./RandNkIext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $Tme $Tmi
	fi

	echo I$I $DIRECTORY "Running ..."
    else
    echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
