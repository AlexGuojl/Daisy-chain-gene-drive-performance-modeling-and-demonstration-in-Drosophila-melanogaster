#!/bin/bash
#SBATCH -J gjl_1d_emb
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_con_sup_%j.out
#SBATCH -e gjl_con_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p emb_res_data

drop_size=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99);

emb_res=(0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);


for ((i=1;i<11;i++)) #repeat times
do
    for ((j=7;j<8;j++)) #j traverse emb_res
    do
        for ((m=11;m<12;m++))  
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            if (($k==111111));then
            python3 daisy_1d_driver_emb.py -drop_size ${drop_size[$m]} -emb_res ${emb_res[$j]} -header > emb_res_data/$k.part
            else
            python3 daisy_1d_driver_emb.py -drop_size ${drop_size[$m]} -emb_res ${emb_res[$j]} > emb_res_data/$k.part
            fi
        done
    done
done


wait
cd emb_res_data


cd ..

