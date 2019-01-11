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
import patsy
import statsmodels.api as sm
import argparse

#function to run the regression for each snp
def doStats(i, genotypes_01, samples):
    n = genotypes_01.shape[1]
    #if allele is present in at least 5 samples
    if (genotypes_01[i,:].sum()) < 5:
        return pd.DataFrame({'Pop_Dev': ['NA'],
                             'Pop_PCs_Dev': ['NA'],
                             'Pop_PCs_Q_Dev': ['NA']})

    elif (genotypes_01[i,:].sum()) > n-5:
        return pd.DataFrame({'Pop_Dev': ['NA'],
                             'Pop_PCs_Dev': ['NA'],
                             'Pop_PCs_Q_Dev': ['NA']})

    try: #subset the genotypes at position i
        samples['genotypes'] =\
                        pd.Series(
                        genotypes_01[i,:],
                        index=samples.index)
        # Population
        y,X_pop = patsy.dmatrices(
                        "genotypes ~ Pop",
                        samples)
        logistic_model_pop =\
                        sm.GLM(y, X_pop,
                        family=sm.families.Binomial())
        results_pop = logistic_model_pop.fit()

        # Population + PCs
        y,X_pop_PCs = patsy.dmatrices(
                            "genotypes ~ Pop+PC1+PC2",
                            samples)
        logistic_model_pop_PCs =\
                            sm.GLM(y, X_pop_PCs,
                            family=sm.families.Binomial())
        results_pop_PCs = logistic_model_pop_PCs.fit()

       # Population + PCs + Q
        y,X_pop_PCs_Q = patsy.dmatrices(
                            "genotypes ~ Pop+PC1+PC2\
                            +average_quality_of_mapped_bases",
                            samples)
        logistic_model_pop_PCs_Q =\
                            sm.GLM(y, X_pop_PCs_Q,
                            family=sm.families.Binomial())
        results_pop_PCs_Q = logistic_model_pop_PCs_Q.fit()

        #return deviance
        return pd.DataFrame({'Pop_Dev': results_pop.deviance,
                             'Pop_PCs_Dev':results_pop_PCs.deviance,
                             'Pop_PCs_Q_Dev':results_pop_PCs_Q.deviance})
    except:
        return pd.DataFrame({'Pop_Dev': ['Error'],
                             'Pop_PCs_Dev': ['Error'],
                             'Pop_PCs_Q_Dev': ['Error']})

def main(args):
    ## VCF file name
    vcf_file_name = 'Chr' + args.chr + '.vcf.recode.vcf.gz'
    vcf_file_path = os.path.join('/Users','luke','genomes','genomes',
                                        '1kGP_NotFiltered',vcf_file_name)
    tabix_file_path = vcf_file_path + ".tbi"
    h5_file_name = 'chr' + args.chr + '.h5'
    h5_file_path = os.path.join('/Users','luke','genomes','genomes',
                                        '1kGP_NotFiltered',h5_file_name)
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
    #this file has Name, Pop, Qual, and the first 5 Global PCs
    samples_fn = os.path.join('/Users','luke','Documents','PCAperPop',
                                        'Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
    samples = pd.DataFrame.from_csv(samples_fn, sep=' ')
    #Verify that individuals are in same order in both files
    assert (np.array(samples.index) == callset['samples'] ).all()
    #reformat the callset
    genotypes = allel.GenotypeChunkedArray(callset['calldata/GT'])
    #Transform hom/het to 0/1
    genotypes_012 = genotypes.to_n_alt(fill=-1)
    genotypes_01 = genotypes_012.astype(bool).astype(int)
    print('Chr ' + args.chr +' Formatting complete')
    #for loop append results
    n=len(genotypes_01)
    results = pd.DataFrame([])
    for i in range(n): 
        out = doStats(i, genotypes_01, samples)
        results = results.append(out, ignore_index=True)

    print('Chr ' + args.chr +' loop complete, writing to file')

    chrPos = pd.DataFrame({'CHROM' : callset['variants']['CHROM'][:n],
                           'POS' : callset['variants']['POS'][:n],
                           'AF' : np.around(callset['variants']['numalt'][:n]\
                               /callset['variants']['NS'][:n]*100,3)})

    output = pd.concat([chrPos, results], axis = 1)

    fileName= args.o + 'Regression_Chr_' + args.chr + '.csv'
    output.to_csv(fileName, mode = 'w', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Logistic Regression Analysis')
    parser.add_argument('-chr', help = 'chromosome')
    parser.add_argument('-o', help= 'output directory')
    args = parser.parse_args()
    main(args)
