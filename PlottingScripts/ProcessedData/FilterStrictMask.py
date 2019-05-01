import multiprocessing as mp
import sys, os
import subprocess

def Bash_cmd(i):
    #Set Path and vcfName
    path = '/Users/luke/genomes/genomes/hg19/phase3/'
    outPath = '/Users/luke/Documents/Legacy_Misc/Strict/'
    vcfName = '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.'
    Bed = '/Users/luke/genomes/BED_MASKS/20141020.strict_mask.whole_genome.bed'

    In = "{0}ALL.chr{1}{2}vcf.gz".format(path,i,vcfName)
    Out = "{0}Chr{1}_strict".format(outPath,i)

    subprocess.call("/usr/local/bin/bcftools-1.6/bcftools filter \
            --regions-file {0} {1} \
            | /usr/local/bin/bcftools-1.6/bcftools view \
            --drop-genotypes --output-type v \
            --output-file {2}".format(Bed, In, Out), shell=True)

if __name__ == '__main__':
    ListOfChrom=list(range(1,23))
    pool_size=7 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    pool.close()
    pool.join()
