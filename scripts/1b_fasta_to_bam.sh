#!/bin/bash

ref="data/ref.fasta"
fastas="data/truth"
# fastas="data/truth data/calls"
Y='\033[0;33m'
W='\033[0m'

for fasta in $fastas; do
    echo -e "\n${Y}[aligning $fasta]${W}"

    # mm2_presets="map-ont asm10 asm5 asm20"
    # for p in $mm2_presets; do
    #     echo -e "${Y}> minimap2 ${p} ${W}"
    #     minimap2 \
    #         -ax $p \
    #         --eqx \
    #         -t `nproc` \
    #         --end-bonus 5000 \
    #         -z 40000,40000 \
    #         --cap-sw-mem 0 \
    #         $ref \
    #         ${fasta}.fasta |
    #     samtools view \
    #         -@ `nproc` \
    #         -F 2304 \
    #         -b |
    #     samtools sort \
    #         -@ `nproc` |
    #     samtools calmd \
    #         -@ `nproc` \
    #         -b \
    #         -Q \
    #         - \
    #         $ref \
    #     > ${fasta}_mm2-${p}.bam 2> ${fasta}_mm2-${p}.log
    #     samtools index ${fasta}_mm2-${p}.bam
    # done

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

    # # lra
    # source ~/miniconda3/etc/profile.d/conda.sh
    # conda activate
    # lra_presets="ONT CONTIG"
    # for p in $lra_presets; do
    #     echo -e "${Y}> lra ${p} ${W}"
    #     lra global -${p} $ref
    #     lra align \
    #         -${p} \
    #         $ref \
    #         ${fasta}.fasta \
    #         -t `nproc` \
    #         -p s > ${fasta}_lra-${p}.sam 

        # |
        # > samtools view \
        #     -@ `nproc` \
        #     -F 2304 \
        #     -b |
        # samtools sort \
        #     -@ `nproc` |
        # samtools calmd \
        #     -@ `nproc` \
        #     -b \
        #     -Q \
        #     - \
        #     $ref \
        # > ${fasta}_lra-${p}.bam 2> ${fasta}_lra-${p}.log

        # samtools index ${fasta}_lra-${p}.bam
    # done
    # conda deactivate

    # npore (use original project, not lib)
    echo -e "${Y}> npore${W}"
    npore_dir="/home/timdunn/npore"
    . $npore_dir/venv3/bin/activate
    python3 $npore_dir/src/realign.py \
        --bam ${fasta}_orig.bam \
        --ref $ref \
        --out_prefix `pwd`/${fasta}_npore \
        --stats_dir $npore_dir/guppy5_stats 2> ${fasta}_npore.log
    deactivate

    samtools view \
        -@ `nproc` \
        -b \
        ${fasta}_npore.sam |
    samtools sort \
        -@ `nproc` |
    samtools calmd \
        -@ `nproc` \
        -b \
        -Q \
        - \
        $ref \
    > ${fasta}_npore.bam
    samtools index ${fasta}_npore.bam

done
    
