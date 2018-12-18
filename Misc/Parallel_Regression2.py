import multiprocessing as mp
import subprocess

def Bash_cmd(i):
    #print statement to keep track of where you are
    print('Working on Chrom : '+str(i))
    #put your bash command below
    subprocess.call('for pop in GBR FIN CHS PUR CDX CLM IBS PEL PJL KHV ACB GWD ESN BEB MSL STU ITU CEU YRI CHB JPT LWK ASW; \
    do echo ${{pop}} {0}; \
    (ls /Users/luke/Documents/FinalRegression/CHR{0}_NOT.Regression_${{pop}}.csv && echo "File Exists") \
    || Rscript /Users/luke/Documents/QualityPaper/Misc/JustStats.R \
    /Users/luke/Documents/Regression/CHR{0}_NOT.Format_${{pop}}.csv \
    /Users/luke/Documents/FinalRegression/CHR{0}_NOT.Regression_${{pop}}.csv \
    ${{pop}};done'.format(i), shell=True)

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
