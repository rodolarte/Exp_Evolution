#PBS -l nodes=1:ppn=1,mem=10gb,walltime=48:00:00 -j oe -N snpEff
#PBS -d ./ -M rolarte@ucr.edu -m abe

module load sickle
module load bwa
module load samtools
module load repeat-masker
module load vcftools
module load picard
module load GATK/3.3.0

JAVA=/opt/java/jdk1.7.0_17/bin/java

list=bam.list

genome=/shared/stajichlab/projects/Hortaea_werneckii/assemblies/Hw2/Hw2.fasta

repeat=$genome.out.bed

prefix=50Hw
species=yeast

mem=1
cpu=1



###Annotation
REFERENCE="Hw2"
DB=/rhome/rolarte/bigdata/Exp_Evolution/snpEff/snpeff
BIN=/rhome/rolarte/bigdata/Exp_Evolution/snpEff/snpeff/snpEff
INPUT_VCF=$prefix.gatk.snp.hardfilter.pass.vcf
VCF_BASE=`basename $INPUT_VCF .vcf`
SNPEFF="java -Xmx1g -jar $BIN/snpEff.jar"
if [ ! -e $DB/$REFERENCE ]; then
echo "Downloading DB $REFERENCE"
$SNPEFF download -c $DB/snpEff.config -v $REFERENCE
fi

if [ ! -e $prefix.gatk.snp.hardfilter.pass.anno.vcf ]; then

$SNPEFF eff -c $BIN/snpEff.config -v -o gatk $REFERENCE $INPUT_VCF > $VCF_BASE.eff.vcf

$JAVA -Xmx10g -jar $GATK \
      -T VariantAnnotator \
      -R $genome \
      -A SnpEff \
      --variant $prefix.gatk.snp.hardfilter.pass.vcf \
      --snpEffFile $prefix.gatk.snp.hardfilter.pass.eff.vcf \
      -o $prefix.gatk.snp.hardfilter.pass.anno.vcf
fi
