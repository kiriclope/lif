#!/bin/bash 

if [ -z "$1" ]; then 
    echo "Network"
    read -p "dir nbpop N K g" dir nbpop N K g
else
    read dir nbpop N K g <<<  "$1 $2 $3 $4 $5"    
    read CrecE CrecI CrecS CrecV Cff <<< "${6} ${7} ${8} ${9} ${10}"
    read IdxMin dIdx IdxMax <<< "${11} ${12} ${13}"
fi

echo IdxMin$IdxMin dIdx$dIdx IdxMax$IdxMax

for Idx in $(seq $IdxMin $dIdx $IdxMax); do
    for I in $(seq .0 .05 .05); do 
	screen -dmS LIF_${dir}_RND_${Idx}_N${N}K${K} ./Jrobust.out ${dir}_RND_$Idx $nbpop $N $K $g $I $CrecE $CrecI $CrecS $CrecV $Cff 
	echo LIF_${dir}_RND_${Idx}_N${N}K${K} "Running ..."
    done
done
