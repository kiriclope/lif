#!/bin/bash
# LANG=en_US

# l_N=(4)
# l_K=(500 1000 2000 3000)

l_N=(8)
l_K=(500 1000)

if [ -z "$1" ]
then 
    echo "Network" 
    read -p "dir nbpop g" dir nbpop g 
else  
    
    read dir nbpop g <<<  "$1 $2 $3" 
    
    if ! [ -z "$4" ]
    then
	if [ "$nbpop" = 2 ]
	then
	    read CrecE CrecI Cff <<<  "$4 $5 $6"
	else
	    read CrecE CrecI CrecS CrecV Cff <<<  "$4 $5 $6 $7 $8"
	fi
    fi
fi

echo Nlist "${l_N[@]}"
echo Klist "${l_K[@]}"

for N in "${l_N[@]}"; do
    for K in "${l_K[@]}"; do
	
	# if [ -z "$4" ]
	# then
	screen -dmS LIF_N${N}K${K} ./ifnet.out $dir $nbpop $N $K $g
	echo LIF_${dir}N${N}K${K} "Running ..."
	# else	    
	#     if [ "$nbpop" = 2 ]
	#     then
	# 	screen -dmS LIF_${dir}N${N}K${K} ./ifnet.out $dir $nbpop $N $K $g $CrecE $CrecI $Cff
	# 	echo LIFSpace_${dir}N${N}K${K}g${g}CrecE${CrecE}CrecI${CrecI}Cff${Cff} "Running ..."
	#     else
	# 	screen -dmS LIF_${dir}N${N}K${K} ./ifnet.out $dir $nbpop $N $K $g $CrecE $CrecI $CrecS $CrecV $Cff	
	# 	echo LIFSpace_${dir}N${N}K${K}g${g}CrecE${CrecE}CrecI${CrecI}CrecS${CrecS}CrecV${CrecV}Cff${Cff} "Running ..."
	#     fi	    
	# fi
	
	# if [ $N=16 ] ; then
	    
	#     screen -dmS LIF_${dir}N${N}K${K} ./large.out $dir $nbpop $N $K $g
	#     echo LIF_${dir}N${N}K${K} "Running ..."
	    
	# else 
	    
	#     screen -dmS LIF_${dir}N${N}K${K} ./a.out $dir $nbpop $N $K $g
	#     echo LIF_${dir}N${N}K${K} "Running ..."
	    
	# fi
	
    done
    
done
