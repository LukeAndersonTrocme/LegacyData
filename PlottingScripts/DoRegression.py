##Run Regression in Python

import allel
import os
import gzip
print(allel.__version__)
import numpy as np
import scipy
import pandas as pd
import h5py
import allel; print('scikit-allel', allel.__version__)
import statsmodels.api as sm
import argparse
import rpy2
print(rpy2.__version__)
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
tuple(ro.globalenv.keys())
import time

#function to run the regression for each snp
def doStats(chunk,samples):
     try:
         ro.globalenv['chunky'] = chunk
         ro.globalenv['samply'] = samples
         #samples['GT'] = np.array(row)

         #formula='GT ~ Pop + PC1 + PC2 + average_quality_of_mapped_bases'
         out = ro.r('apply(chunky, 1, function(x) \
                     anova(glm(x ~ samply$Pop + \
                     samply$PC1 + samply$PC1 + \
                     samply$average_quality_of_mapped_bases,\
                     family=binomial),test="Chi")[5,2])')

         #dev=ro.r('sapply(out, function(x) anova(x, test=Chi)[4,2]')

         #model_pop_PCs_Q = ro.r.glm(formula, data=samples,family='binomial')

         #deviance = ro.r.anova(model_pop_PCs_Q, test='Chi')

         #out = str(deviance[1][4])

         return out
     except Exception as e:
         return str(e)

def main(args):
    #this file has Name, Pop, Qual, and the first 5 Global PCs
    samples_fn = os.path.join('/Users','luke','Documents','PCAperPop',
                                        'Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
    samples = pd.DataFrame.from_csv(samples_fn, sep=' ')

    ## VCF file name
    path = os.path.join('/Users','luke','genomes','genomes','hg19','phase3')
    vcf_file_name = 'ALL.chr' + args.chr +\
            '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    vcf_file_path = os.path.join(path,vcf_file_name)
    tabix_file_path = vcf_file_path + ".tbi"
    h5_file_name = 'chr' + args.chr + '.h5'
    h5_file_path = os.path.join(path,h5_file_name)
    #read h5 file (much faster than vcf)
    exists = os.path.isfile(h5_file_path)
    if exists:
        callset = h5py.File(h5_file_path, mode='r')
        print('Chr ' + args.chr +' H5 exists')
    else:
        callset = allel.read_vcf(vcf_file_path)
        # Uncomment to re-create hdf5 file
        print('Chr ' + args.chr +' H5 being produced')
        allel.vcf_to_hdf5(vcf_file_path, h5_file_path, fields='*', overwrite=True)

    print('Chr ' + args.chr +' File is read')
    start = time.time()
    #Verify that individuals are in same order in both files
    assert (np.array(samples.index) == callset['samples'] ).all()
    genotypes = allel.GenotypeChunkedArray(callset['calldata/GT'])
    #get biallelic sites
    allele_counts = genotypes.count_alleles()
    is_biallelic = allele_counts.is_biallelic_01()

    genotypes_biallelic = genotypes.compress(is_biallelic)
    allele_counts_biallelic_all_alleles = allele_counts.compress(is_biallelic)
    pos_biallelic = callset['variants/POS'][:].compress(is_biallelic)
    #compute allele frequency
    relevant_column = np.array([False] * allele_counts_biallelic_all_alleles.shape[1])
    relevant_column[0:2] = True
    allele_counts_biallelic = allele_counts_biallelic_all_alleles.compress(relevant_column, axis = 1)
    alt_allele_freqs = allele_counts_biallelic[:,1] / allele_counts_biallelic[:].sum(axis = 1)

    #select variants of interest
    variants_selection=(alt_allele_freqs > 3/2504) & (alt_allele_freqs < 2501/2504)
    variants_pass = genotypes_biallelic.compress(variants_selection)
    #Transform hom/het to 0/1
    genotypes_012 = variants_pass.to_n_alt(fill=-1)
    genotypes_01 = genotypes_012.astype(bool).astype(int)
    print('Chr ' + args.chr +' Formatting complete')
    end = time.time()
    print('time for formatting : '+ str(end-start) + 's')
    #for loop append results
    n = 10 #len(genotypes_01)
    s = 100
    iter = 0
    for i in range(n):
        start = time.time()
        chunk = genotypes_01[ iter : iter+s]
        out = doStats(chunk, samples)
        iter = iter+s

        fileName= args.o + 'Chr' + args.chr + '_deviance.csv'
        if not os.path.isfile(fileName):
            f=open(fileName,'w+')
            f.write(str(out) + '\n')
            end = time.time()
            print('time for loop : '+ str(end-start) + 's')
        else:
            f=open(fileName,'a')
            f.write(str(out) + '\n')
            end = time.time()
            print('time for loop : '+ str(end-start) + 's')

    print('Chr ' + args.chr +' loop complete, writing file')

    fileName= args.o + 'Chr' + args.chr + '_deviance.csv'
    results = pd.read_table(fileName,header=None)

    chrPos = pd.DataFrame({'CHROM' : callset['variants']['CHROM'][:].compress(variants_selection),
                           'POS' : callset['variants']['POS'][:].compress(variants_selection),
                           'AF' : alt_allele_freqs.compress(variants_selection)})

    output = pd.concat([chrPos, results], axis = 1)

    fileName= args.o + 'Regression_Chr_' + args.chr + '.csv'
    output.to_csv(fileName, mode = 'w', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Logistic Regression Analysis')
    parser.add_argument('-chr', help = 'chromosome')
    parser.add_argument('-o', help= 'output directory')
    args = parser.parse_args()
    main(args)
