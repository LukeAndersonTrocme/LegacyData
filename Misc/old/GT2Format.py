import multiprocessing as mp
import sys, os
import time

def Bash_cmd(i):
    #print statement to keep track of where you are
    print('Working on Chrom : '+str(i))
    #put your bash command below
    #path='/Volumes/gravel/luke_projects/1000Genomes'
    path='/Users/luke/Documents/Regression/'
    os.system('python ~/Documents/QualityPaper/Misc/GT_wideLong.py \
    -i {0}/CHR{1}.Genotypes_indels.txt.gz \
    -o {0}/CHR{1}_indels.Format'.format(path, i)) #&& \
    #python ./Parallel_Regression.py {1}'.format(path, i))

if __name__ == '__main__':
    #you can replace ListOfChrom with whatever you are iterating through
    ListOfChrom=list(range(1,22))
    #pool_size is number of CPU
    pool_size=6 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    #this is where you are calling the parallelization
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    #some magic, not sure what this does
    pool.close()
    pool.join()
