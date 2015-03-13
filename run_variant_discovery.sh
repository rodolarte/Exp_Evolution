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

#bam=50Hw.RG.bam

prefix=50Hw
species=yeast

mem=1
cpu=1



###snp call
if [ ! -e $prefix.gatk.snp.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
      -T UnifiedGenotyper \
      --filter_mismatching_base_and_quals \
      -R $genome \
      -I $list \
      -o $prefix.gatk.snp.raw.vcf \
      -nct $cpu \
      -stand_call_conf 30 \
      -stand_emit_conf 10 \
      -glm SNP \
      > $prefix.gatk.snp.log 2> $prefix.gatk.snp.log2

fi

###indel call
if [ ! -e $prefix.gatk.indel.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
      -T UnifiedGenotyper \
      --filter_mismatching_base_and_quals \
      -R $genome \
      -I $list \
      -o $prefix.gatk.indel.raw.vcf \
      -nct $cpu \
      -stand_call_conf 30 \
      -stand_emit_conf 10 \
      -glm INDEL \
      > $prefix.gatk.indel.log 2> $prefix.gatk.indel.log2

fi


###hardfilter indel
if [ ! -e $prefix.gatk.indel.hardfilter.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.indel.raw.vcf \
      -o $prefix.gatk.indel.hardfilter.vcf \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "ReadPosRankSum < -20.0" \
      --filterName "ReadPosFilter" \
      --filterExpression "FS > 200.0" \
      --filterName "FSFilter" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "QualFilter"

$JAVA -Xmx1g -jar $GATK \
      -T SelectVariants \
      -R $genome \
      --variant $prefix.gatk.indel.hardfilter.vcf \
      -o $prefix.gatk.indel.hardfilter.pass.vcf \
      --excludeFiltered

fi

###hardfilter snp
if [ ! -e $prefix.gatk.snp.hardfilter.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.raw.vcf \
      -o $prefix.gatk.snp.hardfilter.vcf \
      --clusterSize 3 \
      --clusterWindowSize 10 \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "MQ < 40.0" \
      --filterName "MQFilter" \
      --filterExpression "FS > 60.0" \
      --filterName "FSFilter" \
      --filterExpression "HaplotypeScore > 13.0" \
      --filterName "HaplotypeScoreFilter" \
      --filterExpression "MQRankSum < -12.5" \
      --filterName "MQRankSumFilter" \
      --filterExpression "ReadPosRankSum < -8.0" \
      --filterName "ReadPosRankSumFilter" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "StandardFilters" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --mask $prefix.gatk.indel.hardfilter.pass.vcf \
      --maskName "INDEL"

$JAVA -Xmx1g -jar $GATK \
      -T SelectVariants \
      -R $genome \
      --variant $prefix.gatk.snp.hardfilter.vcf \
      -o $prefix.gatk.snp.hardfilter.pass.vcf \
      --excludeFiltered

fi



