#PBS -q js -l nodes=1:ppn=16 -q highmem -N GATK -j oe -o GATK.log
#PBS -d ./ -M rolarte@ucr.edu -m abe

module load vcftools

##Print summary of snps 

vcf-to-tab < ALL.gatk.snp.hardfilter.pass.eff.vcf > ALL.gatk.snp.hardfilter.pass.eff.tab

##Print stats of snps

if [ -f ALL.gatk.snp.hardfilter.pass.eff.vcf ]; then
 # compress and index
bgzip ALL.gatk.snp.hardfilter.pass.eff.vcf
tabix ALL.gatk.snp.hardfilter.pass.eff.vcf.gz
fi
vcf-stats ALL.gatk.snp.hardfilter.pass.eff.vcf.gz > ALL.gatk.snp.hardfilter.pass.eff.vcfstats




