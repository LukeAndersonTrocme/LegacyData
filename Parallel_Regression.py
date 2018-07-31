import multiprocessing as mp
import sys, os
import time

chrom = sys.argv[1]
global chrom

def Bash_cmd(pop):
    print('Working on Pop : '+str(pop))
    #put your bash command below
    input='~/Documents/Regression/CHR{0}.Format_{1}.csv.gz'.format(chrom,pop)
    output='~/Documents/Regression/CHR{0}.Regression_{1}.csv'.format(chrom,pop)

    print(input)
    os.system('Rscript ~/Documents/QualityPaper/JustStats.R {0} {1} {2}'.format(input, output, pop))

if __name__ == '__main__':
    with open('/Users/luke/genomes/genomes/PopNames/NamePop.txt') as f:
        Pops = f.read().splitlines()
    pool_size=6 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, Pops)
    pool.close()
    pool.join()
