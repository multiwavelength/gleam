#!/bin/bash
#path1='/home/andra/Desktop/Keep/Cluster_spectroscopy/ACRes/line_measurements/allfluxes/'
path1='/media/andra/DATA/Keep/Cluster_spectroscopy/ACRes/line_measurements/allfluxes/'
local=$path1'Z2089/Q6'
local1=$local'/EXT1'
echo $local1

q=0
for i in $(ls $local1/spec1d*); do
    echo $i; 
    s="$(cut -d'.' -f3 <<<$i)"; 
    c="$(cut -d'.' -f2 <<<$i)";
    d="$(cut -d'.' -f4 <<<$i)";
    
    q=$[q+1]
    echo $q
    if [ $q -gt 8 ];
        then
            while pgrep -x "eog" > /dev/null;
                do sleep 1;
            done
            echo 'Large number'
            q=0   
        else
            echo 'Small number'
            eog $local1/linefits.$c.$s.$d.png -n &
    fi
done

