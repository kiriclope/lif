#!/bin/bash 

if [ -z "$1" ]; then 
    echo "Network" 
    read -p "dir nbpop N K g" dir nbpop N K g 
else
    read dir nbpop N K g <<<  "$1 $2 $3 $4 $5" 
    if [ -z "$9" ]; then 
	filename='list.txt'
    else
	read IdxMin dIdx IdxMax Iprtr <<< "${6} ${7} ${8} ${9}"
	echo IdxMin$IdxMin dIdx$dIdx IdxMax$IdxMax
    fi
fi

# if [ -z "$9" ]; then 
#     Iprtr=.05;
#     while read Idx; do 
# 	for I in $(seq 0 $Iprtr $Iprtr); do 
# 	    screen -dmS LIF_${dir}_RND_${Idx}_N${N}K${K}I$I ./Jrobust.out ${dir}/RND/$Idx $nbpop $N $K $g $I 
# 	    echo LIF_${dir}/RND/${Idx}_N${N}K${K}I$I "Running ..."
# 	done
#     done < $filename
# else
    for Idx in $(seq $IdxMin $dIdx $IdxMax); do
	for I in $(seq 0.0 .05 $Iprtr); do 
	    # screen -dmS LIF_${dir}_RND_${Idx}_N${N}K${K}I$I ./Jrobust.out ${dir}/RND/$Idx $nbpop $N $K $g $I 
	    screen -dmS LIF_${dir}_RND_${Idx}_N${N}K${K}I$I ./JrobustS1.out ${dir}/RND/$Idx $nbpop $N $K $g $I 
	    echo LIF_${dir}_RND_${Idx}_N${N}K${K}I$I "Running ..."
	done
    done    
# fi
