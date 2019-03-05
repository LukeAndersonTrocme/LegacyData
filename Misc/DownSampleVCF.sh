##Download 1kGP data and Down Sample to Omni Chip Positions

GenoPath='/Users/luke/genomes/genomes/hg19/phase3'
cd $GenoPath

#GetOmni Chip Positions
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz

gzcat $GenoPath/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz | \
grep -v '^#' | \
awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2}' > \
$GenoPath/ALL.chip.POS

for f in `seq 1 22`
  do echo $f
  #download VCF if needed
  if [ ! -f $GenoPath/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz]; then
    echo "Downloading VCF CHR ${f}"
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  fi
  #Get header of VCF to keep sample names
  if [ $f = 1 ]; then
    grep '^#' $GenoPath/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    > $GenoPath/ALL.chip.vcf
  fi
  #select positions on Omni Chip
  echo "extracting positions for CHR ${f}"
  bcftools view -R $GenoPath/ALL.chip.POS \
    $GenoPath/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
    grep -v '#' >> $GenoPath/ALL.chip.vcf

  done
