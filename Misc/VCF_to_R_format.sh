
for f in `seq 1 21`
do echo $f
java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar \
extractFields $path/ALL.chr${f}. CHROM POS AF "GEN[*].GT" \
| sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d'\
> /Users/luke/genomes/genomes/hg19/Genotypes/CHR${f}.Genotypes.txt; done


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
    extractFields {1}/ALL.chr{2}.{3} CHROM POS AF "GEN[*].GT" \
    | sed "s/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d"\
    | tail -n +2 > /Users/luke/genomes/genomes/hg19/Genotypes/CHR${f}.Genotypes.txt'.format(path,i,name))

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
