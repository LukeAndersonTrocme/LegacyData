import multiprocessing as mp
import sys, os
import time

def Bash_cmd(i):
    #print statement to keep track of where you are
    print('Working on Chrom : '+str(i))
    #put your bash command below
    path='/Users/luke/genomes/genomes/hg19/phase3'
    name='phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

    os.system('java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar \
    extractFields {0}/ALL.chr{1}.{2} CHROM POS AF "GEN[*].GT" \
    | sed "s/0|0/0/g ; s/0|[1-9]/1/g ; s/[1-9]|0/1/g ; s/[1-9]|[1-9]/1/g ; /|/d"\
    | tail -n +2 | awk "\$3 > 0.000599" | gzip > /Users/luke/genomes/genomes/hg19/Genotypes/DoubleCheck/CHR{1}.Genotypes.txt.gz'.format(path,i,name))

if __name__ == '__main__':
    #you can replace ListOfChrom with whatever you are iterating through
    ListOfChrom=list(range(22,0,-1))
    #pool_size is number of CPU
    pool_size=4 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    #this is where you are calling the parallelization
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    #some magic, not sure what this does
    pool.close()
    pool.join()
