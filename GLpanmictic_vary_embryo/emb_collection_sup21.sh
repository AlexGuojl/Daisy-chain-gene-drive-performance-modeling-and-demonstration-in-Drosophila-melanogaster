#!/bin/bash
#SBATCH -J gjl_emb_sup
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o gjl_emb_sup_%j.out
#SBATCH -e gjl_emb_sup_%j.err
#SBATCH --no-requeue
#SBATCH -A jchamper_g1
#SBATCH --qos=jchampercnl
#SBATCH -c 1
pkurun  sleep 1

mkdir -p emb_data_sup

drop_size=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);
emb_res=(0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0);


for ((i=1;i<11;i++)) #go through drop_size 不用改
do
    for ((j=21;j<22;j++)) #j traverse embryo_resistance_rate 每次都需改
    do
        for ((m=1;m<22;m++))  # run i-j 10 times repeat!!,不用改
        do
            k=$((($i+10)*10000+($j+10)*100+$m+10)) # cautious about file order 1 > 02 110000+11*100+11
            python3 daisy_emb_sup_driver_0221.py -drop_size ${drop_size[$m]} -embryo_resistance_rate ${emb_res[$j]} > emb_data_sup/$k.part
            #fi
        done
    done
done


wait
cd emb_data_sup
#cat *.part > emb_nonGLU.csv #每次都要改！
#rm *.part


cd ..