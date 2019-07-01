#!/bin/bash 

if [ -z "$1" ]; then 
    echo "Network"
    read -p "dir nbpop N K g" dir nbpop N K g
else
    read dir nbpop N K g <<<  "$1 $2 $3 $4 $5"    
    read CrecE CrecI CrecS CrecV Cff <<< "${6} ${7} ${8} ${9} ${10}"
    filename="list.txt"
    read IdxMin dIdx IdxMax <<< "${11} ${12} ${13}"
    echo IdxMin$IdxMin dIdx$dIdx IdxMax$IdxMax
fi

# # if [ -z "$11" ]; then
# while read Idx; do
#     for I in $(seq .0 .05 .0); do 
# 	screen -dmS ${dir}_RND_${Idx}_N${N}K${K}I$I ./Jrobust2D.out ${dir}/RND/$Idx $nbpop $N $K $g $I $CrecE $CrecI $CrecS $CrecV $Cff 
# 	echo LIF_${dir}/RND/${Idx}_N${N}K${K}I$I "Running ..." 
#     done
# done < "$filename"
# else
for Idx in $(seq $IdxMin $dIdx $IdxMax); do
    for I in $(seq .05 .05 .05); do 
	screen -dmS ${dir}_RND_${Idx}_N${N}K${K}I$I ./Jrobust2D.out ${dir}/RND/$Idx $nbpop $N $K $g $I $CrecE $CrecI $CrecS $CrecV $Cff 
	echo LIF_${dir}/RND/${Idx}_N${N}K${K}I$I "Running ..."
    done
done
# fi
