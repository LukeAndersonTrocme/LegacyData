import multiprocessing as mp
import subprocess

def Bash_cmd(i):
    #print statement to keep track of where you are
    print('Working on Chrom : '+str(i))
    #put your bash command below
    subprocess.call('python /Users/luke/Documents/QualityPaper/Misc/Make_h5.py -chr {0}'.format(i), shell=True)

if __name__ == '__main__':
    #you can replace ListOfChrom with whatever you are iterating through
    ListOfChrom=list(range(1,22))
    #print(ListOfChrom)
    #pool_size is number of CPU
    pool_size=1   #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    #this is where you are calling the parallelization
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    #some magic, not sure what this does
    pool.close()
    pool.join()
