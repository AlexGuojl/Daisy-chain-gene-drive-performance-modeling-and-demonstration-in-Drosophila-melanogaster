#!/bin/bash
#SBATCH -J gjl_con_sup
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_con_sup_%j.out
#SBATCH -e gjl_con_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p con_data_sup

drop_size=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);
con_rate=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);


for ((i=1;i<11;i++)) #repeat times
do
    for ((j=16;j<17;j++)) #j traverse embryo_resistance_rate 每次都需改
    do
        for ((m=1;m<22;m++))  #  traverse drop_size 
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            if (($k==111111));then
            python3 daisy_con_sup_driver_0229.py -drop_size ${drop_size[$m]} -conversion_rate ${con_rate[$j]} -header > con_data_sup/$k.part
            else
            python3 daisy_con_sup_driver_0229.py -drop_size ${drop_size[$m]} -conversion_rate ${con_rate[$j]} > con_data_sup/$k.part
            fi
        done
    done
done


wait
cd con_data_sup

cd ..