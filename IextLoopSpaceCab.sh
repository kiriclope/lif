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
	read E I S V  Tme Tmi Tms Tmv Crec Cff dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21}"
    fi

    if [ "$nbpop" = 3 ]; then 
	read E I S <<< "${10} ${11} ${12}"
	read Tme Tmi Tms <<< "${13} ${14} ${15}"
	read CrecEE CrecEI CrecES <<<  "${16} ${17} ${18}"
	read CrecIE CrecII CrecIS <<< "${19} ${20} ${21}"
	read CrecSE CrecSI CrecSS <<< "${22} ${23} ${24}"
	read Cff dI Imax <<< "${25} ${26} ${27}"
    fi

    if [ "$nbpop" = 2 ]; then 
	read E I Tme Tmi CrecEE CrecEI CrecIE CrecII Cff dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20}"
    fi

    if [ "$nbpop" = 1 ]; then 
	read I Tmi CrecI Cff dI Imax <<< "${10} ${11} ${12} ${13} ${14} ${15}"
    fi

fi

export $LANG=en_US

Imin=$I
echo Imin$Imin dI$dI Imax$Imax
for I in $(seq $Imin $dI $Imax); do
# for I in .5 5 10 25 50 75 100 125 150 200 ; do
    
    DIRECTORY=$(printf "../%dpop/Network/%s/N%d/K%d/g%.2f/Prtr_I/Iext_I%.4f/Space/CrecE%.2fCrecI%.2fCff%.2f" "$nbpop" "$dir" "$N" "$k" "$g" "$I" "$CrecE" "$CrecI" "$Cff")
    
    if [ ! -d "$DIRECTORY" ]; then
        
	if [ "$nbpop" = 4 ]; then 
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./ifspace_Iext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $S $V $Tme $Tmi $Tms $Tm $Crec $Cff
	fi
	
	if [ "$nbpop" = 3 ]; then 
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./ifspaceCab.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $S $Tme $Tmi $Tms $CrecEE $CrecEI $CrecES $CrecIE $CrecII $CrecIS $CrecSE $CrecSI $CrecSS $Cff
	fi
	
	if [ "$nbpop" = 2 ]; then 
	    # screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./ifspace_Iext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $Tme $Tmi $CrecE $CrecI $Cff
	    screen -dmS IF${nbpop}pop${dir}N${N}K${k}g${g}I${I} ./GaussCabNkIext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $E $I $Tme $Tmi $CrecEE $CrecEI $CrecIE $CrecII $Cff
	fi
	
	if [ "$nbpop" = 1 ]; then 
	    screen -dmS ${nbpop}pop_${dir}_K${k}I${I}CrecI${CrecI}Cff${Cff} ./ifspace_Iext.out $dir $nbpop $N $k $dt $duration $Tst $Tw $g $I $Tmi $CrecI $Cff
	fi
	
	echo ${nbpop}pop_${dir}_K${k}I${I}CrecI${CrecI}Cff${Cff} $DIRECTORY "Running ..."
    else
	echo I$I $DIRECTORY "Folder already exists ..."
    fi

done
