import multiprocessing as mp
import sys, os

def Bash_cmd(i):
    #Set Path and vcfName
    path='/Users/luke/genomes/genomes/hg19/'

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
