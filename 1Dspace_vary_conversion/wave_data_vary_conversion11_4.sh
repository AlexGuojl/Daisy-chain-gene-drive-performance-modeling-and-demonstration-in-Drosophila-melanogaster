#!/bin/bash
#SBATCH -J gjl_1d_dc
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_con_sup_%j.out
#SBATCH -e gjl_con_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p conversion_vs_dropsize

drop_size=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99);

drive_conversion=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);


for ((i=1;i<11;i++)) #repeat times
do
    for ((j=4;j<5;j++)) #conversion
    do
        for ((m=11;m<12;m++))  #
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            if (($k==111111));then
            python3 daisy_1d_driver_con.py -drop_size ${drop_size[$m]} -conversion_rate ${drive_conversion[$j]} -header > conversion_vs_dropsize/$k.part
            else
            python3 daisy_1d_driver_con.py -drop_size ${drop_size[$m]} -conversion_rate ${drive_conversion[$j]} > conversion_vs_dropsize/$k.part
            fi
        done
    done
done


wait
cd conversion_vs_dropsize


cd ..

