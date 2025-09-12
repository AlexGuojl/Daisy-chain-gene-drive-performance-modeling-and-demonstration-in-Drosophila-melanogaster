#!/bin/bash
#SBATCH -J gjl_som_sup
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_som_sup_%j.out
#SBATCH -e gjl_som_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p som_data_sup

drop_size=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);
som_rate=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);



for ((i=1;i<11;i++)) #repeat times
do
    for ((j=20;j<21;j++)) #j traverse embryo_resistance_rate 每次都需改
    do
        for ((m=1;m<22;m++))  #  traverse drop_size 
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            if (($k==111111));then
            python3 daisy_som_sup_driver_0229.py -drop_size ${drop_size[$m]} -somatic_mult_f ${som_rate[$j]} -header > som_data_sup/$k.part
            else
            python3 daisy_som_sup_driver_0229.py -drop_size ${drop_size[$m]} -somatic_mult_f ${som_rate[$j]} > som_data_sup/$k.part
            fi
        done
    done
done


wait
cd som_data_sup

cd ..