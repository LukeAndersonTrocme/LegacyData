import pandas as pd
import argparse

def main(args):
    #read input file (VCF file cut first 4 cols)
    vcfPos=pd.read_table(args.input,
                         names=["Chrom","Pos", "Ref", "Alt"])

    vcfPos.drop_duplicates(subset='Pos', keep = False)
    #set default Derived Allele to 1
    vcfPos['DerivedAllele'] = 1

    #read chimp human alignment
    chimp=pd.read_table('/Users/luke/bin/smaller_mut_spectrum_pipeline/'
    +'hg19_chimp_align/human_chimp_diffs_chr'
    + args.chrom + '.txt', sep=" ")
    chimp= chimp[chimp['SNP/Indel'] == 'SNP']

    #read oneLine fasta ref file
    ref=open('/Users/luke/bin/smaller_mut_spectrum_pipeline/'
    +'hg19_reference/chr' + args.chrom + '_oneline.txt')
    refseq=ref.read()

    #get context (i.e. one base up- and downstream)
    vcfPos['Context'] = vcfPos.apply(\
    lambda row : refseq[row.Pos - 2 : row.Pos + 1], axis=1)
    vcfPos=vcfPos.drop_duplicates(subset='Pos')

    #only keep SNPs
    DNA=['A','C','G','T']
    vcfPos= vcfPos[vcfPos['Ref'].isin(DNA)]
    vcfPos= vcfPos[vcfPos['Alt'].isin(DNA)]
    vcfPos= vcfPos[~vcfPos['Context'].isin(['N'])]

    #get sites shared by both files
    iPos = set(chimp['Pos']) & set(vcfPos['Pos'])
    iChimp = chimp[chimp['Pos'].isin(list(iPos))]
    iVCF = vcfPos[vcfPos['Pos'].isin(list(iPos))]
    print('Chrom : '+args.chrom)
    print(len(iVCF), len(iChimp))
    assert len(iVCF) == len(iChimp)
    #find sites where Chimp is same as Human Alt
    ChimpDiff = iVCF[iVCF['Alt'] == list(iChimp['Chimp'])]['Pos']
    #ChimpDiff = vcfPos['Pos'][vcfPos['Pos'].isin(chimp['Pos'])][vcfPos['Alt'] == list(chimp['Chimp'])]['Pos']
    #print(ChimpDiff.head())

    ##deal with an error that happens for some chromosomes
    #ValueError: Arrays were different lengths: 321999 vs 321997

    #For sites that have ChimpDiff, change DerivedAllele, Context, Ref and Alt
    vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), \
    'DerivedAllele'] = 0
    vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), \
    'Context'] = vcfPos['Context'].str[0] \
                + vcfPos['Alt'].str[0] \
                + vcfPos['Context'].str[2]

    vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), \
    'Alt'] = vcfPos['Ref'].str[0]
    vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), \
    'Ref'] = vcfPos['Context'].str[1]

    #write table
    vcfPos.to_csv(args.out + 'Chr' + args.chrom + \
    '.AncestralContext.txt', index = False, header = False, sep="\t")

parser = argparse.ArgumentParser(description = 'get Ancestral Context')
parser.add_argument('-input', help = 'input file')
parser.add_argument('-chrom', help = 'chromosome')
parser.add_argument('-out', help= 'output directory')
args = parser.parse_args()
main(args)
