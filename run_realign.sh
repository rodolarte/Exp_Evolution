#PBS -l nodes=1:ppn=1,mem=10gb,walltime=48:00:00 -j oe -N bwa
#PBS -d ./ -M rolarte@ucr.edu -m abe

module load sickle
module load bwa
module load samtools
module load repeat-masker
module load vcftools
module load picard
module load GATK/3.3.0
#PICARD=/opt/picard/1.81
#GATK=/opt/GATK/3.3.0/GenomeAnalysisTK.jar
JAVA=/opt/java/jdk1.7.0_17/bin/java

list=bam.list

genome=/shared/stajichlab/projects/Hortaea_werneckii/assemblies/Hw2/Hw2.fasta

repeat=$genome.out.bed

bam=50Hw.RG.bam

prefix=50Hw
species=yeast

mem=1
cpu=1


## Creating the fasta sequence dictionary file.  We use CreateSequenceDictionary.jar
## from Picard to create a .dict file from a fasta file


fabase=`basename $genome .fasta`

if [ ! -e $fabase.dict ]; then

$JAVA -Xmx3g -jar $PICARD/CreateSequenceDictionary.jar R=$genome O=$fabase.dict
samtools faidx $genome

fi

## Identify intervals around variants
if [ ! -e $prefix.gatk.realign.bam ]; then
$JAVA -Xmx10g -jar $GATK \
       -T RealignerTargetCreator \
       -R $genome \
       -o $prefix.gatk.intervals \
       -I $bam

## Realign based on these intervals
$JAVA -Xmx10g -jar $GATK \
       -T IndelRealigner \
       -R $genome \
       -targetIntervals $prefix.gatk.intervals \
       -I $bam \
       -o $prefix.gatk.realign.bam

samtools index $prefix.gatk.realign.bam
echo $prefix.gatk.realign.bam > bam.list

fi


