#!/bin/bash

if [ -z "$1" ]; then 
    echo "Network" 
    read -p "dir nbpop N k g" dir nbpop N k g 
    read -p "dI JabMin dJab JabMax" <<< "${6} ${7} ${8} ${9}"
else
    read dir nbpop N k g <<<  "$1 $2 $3 $4 $5"    
    read dI JabMin dJab JabMax <<< "${6} ${7} ${8} ${9}"
fi

echo dI$dI JabMin$JabMin dJab$dJab JabMax$JabMax
for Jab in $(seq $JabMin $dJab $JabMax); do 
    screen -dmS LIF_${dir}_N${N}K${k}Iext${I}Jab${Jab} ./JabLoop.out $dir $nbpop $N $k $g $dI $Jab
    echo "Running " LIF_${dir}_N${N}K${k}Iext${I}Jab${Jab} 
done
