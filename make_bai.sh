#PBS -l nodes=1:ppn=1,mem=1gb -j oe -N bamindex
module load samtools

samtools index E100Hw.RG.bam
