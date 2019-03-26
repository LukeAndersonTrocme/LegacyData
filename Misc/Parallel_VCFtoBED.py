import multiprocessing as mp
import sys, os

def Bash_cmd(i):
    #Set Path and vcfName
    path='/Users/luke/genomes/genomes/hg19/'
    vcfName='.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.'

    #Filter Allele Frequency
    fname = '{0}/AF_VCF/chr{1}_AF.vcf.gz'.format(path, i)
    if os.path.isfile(fname) :
        print('VCF already Filtered : '+str(i))
    else :
        print('AF Filter Chrom : '+str(i))
        os.system('/usr/local/bin/bcftools-1.6/bcftools filter \
        {0}/phase3/ALL.chr{1}{2}vcf.gz \
        -e "MAF > 0.9988019169 || MAF < 0.00119808" \
        -Oz -o {0}/AF_VCF/chr{1}_AF.vcf.gz'.format(path, i,vcfName))

    #Split VCF into 100,000 SNP chunks
    fname = '{0}/AF_VCF/Split_chr{1}_aa'.format(path, i)
    if os.path.isfile(fname) :
        print('VCF already Split : '+str(i))
    else :
        print('Split Chrom : '+str(i))
        os.system('gzcat {0}/AF_VCF/chr{1}_AF.vcf.gz \
        | head -1000 | grep "^#" \
        > {0}/AF_VCF/chr{1}_header ; \
        zgrep -v "^#" {0}/AF_VCF/chr{1}_AF.vcf.gz \
        > {0}/AF_VCF/chr{1}_variants ; \
        split -l 100000 \
        {0}/AF_VCF/chr{1}_variants \
        {0}/AF_VCF/Split_chr{1}_'.format(path, i))

    fname = '{0}/AF_VCF/Split_chr{1}_aa.vcf'.format(path, i)
    if os.path.isfile(fname) :
        print('VCF already Combined : '+str(i))
    else :
        print('Combine Chrom : '+str(i))
        os.system('for i in {0}/AF_VCF/Split_chr{1}_*; \
        do cat {0}/AF_VCF/chr{1}_header $i \
        > $i.vcf && rm -f $i ; done ; \
        rm -f {0}/AF_VCF/chr{1}_header {0}/AF_VCF/chr{1}_variants'.format(path, i))

    #Convert VCF to Plink
    print('Plink Chrom : '+str(i))
    os.system('for j in {0}/AF_VCF/Split_chr{1}_*.vcf; \
    do /Users/luke/bin/plink_mac/plink \
    --vcf $j --make-bed --out ${{j%.vcf}} ; \
    done ; rm -f $j'.format(path, i))

if __name__ == '__main__':
    ListOfChrom=list(range(1,23))
    pool_size=7 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    pool.close()
    pool.join()
