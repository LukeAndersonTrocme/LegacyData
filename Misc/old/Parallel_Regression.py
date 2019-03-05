import multiprocessing as mp
import sys, os
import time

chrom = sys.argv[1]
path = sys.argv[2]
def Bash_cmd(pop):
    print('Working on Pop : '+str(pop))
    #put your bash command below
    #input='{2}/Format/CHR{0}.Format_{1}.csv.gz'.format(chrom, pop, path)
    #output='{2}/Regression/CHR{0}.Regression_{1}.csv'.format(chrom, pop, path)

    input='{2}/CHR{0}_indels.Format_{1}.csv.gz'.format(chrom, pop, path)
    output='{2}/CHR{0}_indels.Regression_{1}.csv'.format(chrom, pop, path)

    print('Rscript ~/Documents/QualityPaper/Misc/JustStats.R {0} {1} {2}'.format(input, output, pop))
    os.system('Rscript ~/Documents/QualityPaper/Misc/JustStats.R {0} {1} {2}'.format(input, output, pop))

if __name__ == '__main__':
    with open('/Users/luke/genomes/genomes/PopNames/NamePop.txt') as f:
        Pops = f.read().splitlines()
    pool_size=6 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, Pops)
    pool.close()
    pool.join()
