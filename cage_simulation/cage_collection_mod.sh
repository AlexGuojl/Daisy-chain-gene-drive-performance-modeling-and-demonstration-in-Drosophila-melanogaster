#!/bin/bash
#SBATCH -J gjl_cage_mode
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_cage_mod_%j.out
#SBATCH -e gjl_cage_mod_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p cage_data_mod

drop_size=(0 0.45);

for ((i=1;i<101;i++)) #repeat times
do
    for ((m=1;m<2;m++))  #  traverse drop_size 
    do
        k=$((($i+10)*10000+($m+10)*100)) # cautious about file order 1 > 02 110000+11*100+11
        if (($k==111100));then
        python3 cage_driver_mod.py -drop_size ${drop_size[$m]} -header > cage_data_mod/$k.part
        else
        python3 cage_driver_mod.py -drop_size ${drop_size[$m]} > cage_data_mod/$k.part
        fi
    done
done


wait
cd cage_data_mod

cd ..