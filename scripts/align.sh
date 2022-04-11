#!/bin/bash

ref="input/ref.fasta"
fastas="input/truth input/calls"
G='\033[0;32m'
Y='\033[0;33m'
W='\033[0m'

for fasta in $fastas; do
    echo -e "\n${G}[aligning $fasta]${W}"

    # mm2-ont
    echo -e "${Y}> mm2-ont${W}"
    minimap2 \
        -ax map-ont \
        --eqx \
        -t `nproc` \
        $ref \
        ${fasta}.fasta |
    samtools view \
        -@ `nproc` \
        -b |
    samtools sort \
        -@ `nproc` |
    samtools calmd \
        -@ `nproc` \
        -b \
        -Q \
        - \
        $ref \
    > ${fasta}_mm2-ont.bam 2> ${fasta}_mm2-ont.log
    samtools index ${fasta}_mm2-ont.bam

    # mm2-asm10
    echo -e "${Y}> mm2-asm10${W}"
    minimap2 \
        -ax map-asm10 \
        --eqx \
        -t `nproc` \
        $ref \
        ${fasta}.fasta |
    samtools view \
        -@ `nproc` \
        -b |
    samtools sort \
        -@ `nproc` |
    samtools calmd \
        -@ `nproc` \
        -b \
        -Q \
        - \
        $ref \
    > ${fasta}_mm2-asm10.bam 2> ${fasta}_mm2-asm10.log
    samtools index ${fasta}_mm2-asm10.bam

    # # ngmlr
    # echo -e "${Y}> ngmlr${W}"
    # source ~/miniconda3/etc/profile.d/conda.sh
    # conda activate
    # ngmlr \
    #     -t `nproc` \
    #     -r $ref \
    #     -q ${fasta}.fasta \
    #     -o ${fasta}_ngmlr.bam 2> ${fasta}_ngmlr.log
    # conda deactivate
    # samtools index ${fasta}_ngmlr.bam

    # npore
    echo -e "${Y}> npore${W}"
    npore_dir="/home/timdunn/npore"
    . $npore_dir/venv3/bin/activate
    python3 $npore_dir/src/realign.py \
        --bam ${fasta}_orig.bam \
        --ref $ref \
        --out_prefix ${fasta}_npore \
        --stats_dir $npore_dir/guppy5_stats 2> ${fasta}_npore.log
    deactivate
    samtools index ${fasta}_npore.bam

done
    
