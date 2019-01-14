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
import rpy2
print(rpy2.__version__)
import rpy2.robjects as robjects
import argparse

#function to run the regression for each snp
def doStats(i, genotypes_01, samples, X_pop_PCs, X_pop_PCs_Q):
    try: #subset the genotypes at position i
        GT = pd.Series(genotypes_01[i,:],index=samples.index)
        # Population + PCs
        logistic_model_pop_PCs =\
                            sm.GLM(GT, X_pop_PCs,
                            family=sm.families.Binomial())
        results_pop_PCs = logistic_model_pop_PCs.fit_regularized(alpha= 1.e-7, L1_wt=0)

        # Population + PCs + Q
        logistic_model_pop_PCs_Q =\
                            sm.GLM(GT, X_pop_PCs_Q,
                            family=sm.families.Binomial())
        results_pop_PCs_Q = logistic_model_pop_PCs_Q.fit_regularized(alpha= 1.e-7, L1_wt=0)

        #return deviance
        return pd.DataFrame({'Pop_PCs_Dev':results_pop_PCs.deviance,
                             'Pop_PCs_Q_Dev':results_pop_PCs_Q.deviance})
    except Exception as e:
        return pd.DataFrame({'Pop_PCs_Dev': [e],
                             'Pop_PCs_Q_Dev': [e]})

def main(args):
    #this file has Name, Pop, Qual, and the first 5 Global PCs
    samples_fn = os.path.join('/Users','luke','Documents','PCAperPop',
                                        'Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
    samples = pd.DataFrame.from_csv(samples_fn, sep=' ')
    #Make Design Matrix
    X_pop_PCs = patsy.dmatrix("Pop + PC1 + PC2",samples)

    X_pop_PCs_Q = patsy.dmatrix("Pop + PC1 + PC2 +\
                                average_quality_of_mapped_bases",samples)
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
    #for loop append results
    n=10000
    results = pd.DataFrame([])
    for i in range(n):
        out = doStats(i, genotypes_01, samples, X_pop_PCs, X_pop_PCs_Q)
        results = results.append(out, ignore_index=True)

    print(results)

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
