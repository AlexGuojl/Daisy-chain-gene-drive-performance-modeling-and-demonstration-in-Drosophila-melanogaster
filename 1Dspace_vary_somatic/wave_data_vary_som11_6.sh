#!/bin/bash
#SBATCH -J gjl_1d_dsom
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_som_sup_%j.out
#SBATCH -e gjl_som_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p size_som

drop_size=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99);
som_mult=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);



for ((i=1;i<11;i++)) #repeat times
do
    for ((j=6;j<7;j++)) #j traverse drop_size 每次都需改
    do
        for ((m=11;m<12;m++))  #  traverse emb_res
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            if (($k==111111));then
            python3 daisy_1d_driver_som.py -drop_size ${drop_size[$m]} -somatic_fitness_f ${som_mult[$j]} -header > size_som/$k.part
            else
            python3 daisy_1d_driver_som.py -drop_size ${drop_size[$m]} -somatic_fitness_f ${som_mult[$j]} > size_som/$k.part
            fi
        done
    done
done


wait
cd size_som
#cat *.part > one_dim_drop_vary_size_radius.csv #每次都要改！
#rm *.part

cd ..

