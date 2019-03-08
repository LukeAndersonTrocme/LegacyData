import multiprocessing as mp
import sys, os
import time


def Bash_cmd(i):
    GenoName='phase3_shapeit2_mvncall_integrated_v5a.20130502'
    #print statement to keep track of where you are
    print('Working on Chrom : '+str(i))
    #put your bash command below
    path='/Users/luke/genomes/genomes/hg19/'
    os.system('/Users/luke/bin/plink_mac/plink \
    --vcf {0}/phase3//ALL.chr{1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --make-bed \
    --out {0}/plink/chr{1}'.format(path, i))

if __name__ == '__main__':
    #you can replace ListOfChrom with whatever you are iterating through
    ListOfChrom=list(range(1,23))
    #pool_size is number of CPU
    pool_size=7 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    #this is where you are calling the parallelization
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    #some magic, not sure what this does
    pool.close()
    pool.join()
