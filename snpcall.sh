#PBS -q js -l nodes=1:ppn=16 -q highmem -N GATK -j oe -o GATK.log
#PBS -d ./ -M rolarte@ucr.edu -m abe


module load sickle
module load bwa
module load samtools
module load repeat-masker
PICARD=/opt/picard/1.81
GATK=/opt/GATK/3.1-1/GenomeAnalysisTK.jar
JAVA=/opt/java/jdk1.7.0_17/bin/java
list=bam.list

#fq1=W303_chrII_1.fastq
fq1=1_TAAGGCGA-TAGATCGC_L002_R1_001.fastq
#fq2=W303_chrII_2.fastq
fq2=1_TAAGGCGA-TAGATCGC_L002_R2_001.fastq
#RGID=W303
RGID=0Hw
RGLB=SRR527545

#genome=`pwd`/Saccharomyces.fa
genome=/shared/stajichlab/projects/Hortaea_werneckii/assemblies/Hw2/Hw2.fasta

repeat=$genome.out.bed

#bam=W303.readgroup.bam
bam=0Hw.readgroup.bam

prefix=0Hw
species=yeast
temp=`pwd`/temp
mem=1
cpu=1


if [ ! -e $temp ]; then

mkdir $temp

fi


## Trim sequences

fq1base=`basename $fq1 .fastq`
fq2base=`basename $fq2 .fastq`


if [ ! -e $RGID.readgroup.bam ]; then

sickle pe -f $fq1 -r $fq2 -t sanger -o $fq1base.trimmed.fastq -p $fq2base.trimmed.fastq -s singles_file.trimmed.fastq

## Index genome before we can align (only need to do this once)
bwa index $genome

bwa aln -q 20 -t 16 -f $fq1base.trimmed.sai $genome $fq1base.trimmed.fastq

bwa aln -q 20 -t 16 -f $fq2base.trimmed.sai $genome $fq2base.trimmed.fastq

bwa sampe -f $RGID.sam $genome $fq1base.trimmed.sai $fq2base.trimmed.sai $fq1base.trimmed.fastq $fq2base.trimmed.fastq

samtools view -b -S $RGID.sam > $RGID.unsrt.bam

samtools sort $RGID.unsrt.bam $RGID.sorted

samtools index $RGID.sorted.bam

## Need to Deduplicate reads
$JAVA -Xmx3g -Djava.io.tmpdir=$temp -jar $PICARD/MarkDuplicates.jar INPUT=$RGID.sorted.bam OUTPUT=$RGID.dedup.bam METRICS_FILE=$RGID.dedup.metrics CREATE_INDEX=true VALIDATION_STRIN
GENCY=SILENT

## Fixing Read-Groups
$JAVA -Xmx3g -Djava.io.tmpdir=$temp -jar $PICARD/AddOrReplaceReadGroups.jar INPUT=$RGID.dedup.bam OUTPUT=$RGID.readgroup.bam SORT_ORDER=coordinate CREATE_INDEX=true RGID=$RGID RGLB=
$RGLB RGPL=Illumina RGPU=Genomic RGSM=$RGID VALIDATION_STRINGENCY=SILENT


fi

echo "ReduceReads for each bam using nway co-reduce, which could keep these reads that not have SNP in some sample in common SNPs site"
$JAVA -Xmx10g -jar $GATK -T ReduceReads -R $genome -I $BAM -o reduced.bam

#echo "Reduce Done"
echo "Multisample call and filter by hardfilter or VQSR"


##
# Creating the fasta sequence dictionary file.  We use CreateSequenceDictionary.jar
# from Picard to create a .dict file from a fasta file

#fabase=`basename $genome .fa`
fabase=`basename $genome .fasta`

if [ ! -e $fabase.dict ]; then

$JAVA -Djava.io.tmpdir=$temp -jar $PICARD/CreateSequenceDictionary.jar R=$genome O=$fabase.dict
samtools faidx $genome

fi

## Identify intervals around variants
if [ ! -e $prefix.gatk.realign.bam ]; then
$JAVA -Xmx10g -Djava.io.tmpdir=$temp -jar $GATK \
       -T RealignerTargetCreator \
       -R $genome \
       -o $prefix.gatk.intervals \
       -I $bam

## Realign based on these intervals
$JAVA -Xmx10g -Djava.io.tmpdir=$temp -jar $GATK \
       -T IndelRealigner \
       -R $genome \
       -targetIntervals $prefix.gatk.intervals \
       -I $bam \
       -o $prefix.gatk.realign.bam
samtools index $prefix.gatk.realign.bam
echo $prefix.gatk.realign.bam > bam.list

fi


###snp call
if [ ! -e $prefix.gatk.snp.raw.vcf ]; then
$JAVA -Xmx10g -Djava.io.tmpdir=$temp -jar $GATK \
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
$JAVA -Xmx10g -Djava.io.tmpdir=$temp -jar $GATK \
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

###repeat bed

if [ ! -e $genome.out.bed ]; then

RepeatMasker -qq --species $species $genome
perl repeat2gff.pl $genome.out > $genome.out.gff
perl orderchr.pl --vcf $prefix.gatk.snp.raw.vcf --gff $genome.out.gff

fi

###hardfilter indel
if [ ! -e $prefix.gatk.indel.pass.vcf ]; then
$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
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

$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
      -T SelectVariants \
      -R $genome \
      --variant $prefix.gatk.indel.hardfilter.vcf \
      -o $prefix.gatk.indel.pass.vcf \
      --excludeFiltered

fi

###hardfilter snp
if [ ! -e $prefix.gatk.snp.pass.vcf ]; then
$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
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
      --mask $prefix.gatk.indel.pass.vcf \
      --maskName "INDEL"

$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
      -T SelectVariants \
      -R $genome \
      --variant $prefix.gatk.snp.hardfilter.vcf \
      -o $prefix.gatk.snp.pass.vcf \
      --excludeFiltered

fi


####mask repeat
#if [ ! -e $prefix.gatk.snp.pass.MASK.pass.vcf ]; then
#$JAVA -Xmx1g -jar $GATK \
#      -T VariantFiltration \
#      -R $genome \
#      --variant $prefix.gatk.snp.pass.vcf \
#      -o $prefix.gatk.snp.pass.MASK.vcf \
#      --mask $repeat \
#      --maskName "REPEAT"
#
#$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
#      -T SelectVariants \
#      -R $genome \
#      --variant $prefix.gatk.snp.pass.MASK.vcf \
#      -o $prefix.gatk.snp.pass.MASK.pass.vcf \
#      --excludeFiltered
#fi
#
####Annotation
#REFERENCE="EF4.70"
#DB=/rhome/rolarte/bigdata/CSHL_NGS/analysis1_RO_10062014/snpeff
#BIN=/rhome/rolarte/bigdata/CSHL_NGS/analysis1_RO_10062014/snpeff/snpEff
#INPUT_VCF=$prefix.gatk.snp.pass.MASK.pass.vcf
#VCF_BASE=`basename $INPUT_VCF .vcf`
#SNPEFF="java -Xmx1g -jar $BIN/snpEff.jar"
#if [ ! -e $DB/$REFERENCE ]; then
#echo "Downloading DB $REFERENCE"
#$SNPEFF download -c $DB/snpEff.config -v $REFERENCE
#fi
#
#if [ ! -e $prefix.gatk.snp.pass.MASK.pass.anno.vcf ]; then
#
#$SNPEFF eff -c $BIN/snpEff.config -v -o gatk $REFERENCE $INPUT_VCF > $VCF_BASE.eff.vcf
#
#$JAVA -Xmx1g -Djava.io.tmpdir=$temp -jar $GATK \
#      -T VariantAnnotator \
#      -R $genome \
#      -A SnpEff \
#      --variant $prefix.gatk.snp.pass.MASK.pass.vcf \
#      --snpEffFile $prefix.gatk.snp.pass.MASK.pass.eff.vcf \
#      -o $prefix.gatk.snp.pass.MASK.pass.anno.vcf
#fi

##Print summary of snps and indels

vcf-to-tab < $prefix.gatk.snp.hardfilter.vcf > $prefix.gatk.snp.hardfilter.tab

vcf-to-tab < $prefix.gatk.indel.hardfilter.vcf > $prefix.gatk.indel.hardfilter.tab
