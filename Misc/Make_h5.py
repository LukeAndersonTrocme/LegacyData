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

def main(args):
    ## VCF file name
    path = os.path.join('/Users','luke','genomes','genomes','hg19','phase3')
    vcf_file_name = 'ALL.chr' + args.chr +\
            '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    vcf_file_path = os.path.join(path,vcf_file_name)
    tabix_file_path = vcf_file_path + ".tbi"
    h5_file_name = 'chr' + args.chr + '.h5'
    h5_file_path = os.path.join(path,h5_file_name)
    print(vcf_file_path)
    print(h5_file_path)

    print('Chr ' + args.chr +' VCF being read')
    #read h5 file (much faster than vcf)
    callset = allel.read_vcf(vcf_file_path)
    print('Chr ' + args.chr +' VCF read')
    # Uncomment to re-create hdf5 file

    allel.vcf_to_hdf5(vcf_file_path, h5_file_path, fields='*', overwrite=True)
    print('Chr ' + args.chr +' H5 produced')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'make h5')
    parser.add_argument('-chr', help = 'chromosome')
    args = parser.parse_args()
    main(args)
