import argparse
import pandas as pd
import os

def main(args):
    
    #Read Quality file
    Qual=pd.read_table('~/Dropbox/LukeTemp/SubMeta.txt', sep=' ', quotechar='"')
    #Read Genotype file
    Chunk = pd.read_table(args.i, header=None, 
                     names=Qual.Name, chunksize=100000, compression = 'gzip')
    for GT in Chunk:
        GT = GT.reset_index().rename(columns = {'level_0':'Chr','level_1':'Pos'})
        print('table read')
        chrPos=pd.Series(['Chr','Pos'])
        #run per Pop to save memory
        for p in Qual.Pop.unique():
            print(p, end=" ", flush=True)
            #ID names per Pop
            names=chrPos.append(Qual.loc[Qual['Pop']==p].Name)
            subGT=GT[names]

            #melt from wide to long
            m=pd.melt(subGT,id_vars=['Chr','Pos'], 
                      value_vars=Qual[Qual.Pop==p].Name).rename(
                        columns={'variable':'Name'})
            #replace HomAlt with Het !!!!
            m.value=m.value.replace(2,1)
            #merge Quality with GT
            GT_Qual = pd.merge(m, Qual[Qual.Pop==p], on='Name')
            del m
            #get number of unique GT per Pop per Pos
            r=GT_Qual.groupby(['Pos'])['value'].nunique().reset_index().rename(columns={'value':'fixed'})

            GT_Qual=pd.merge(GT_Qual, r, on=['Pos'], how='outer')
            #get rid of fixed sites
            GT_Qual=GT_Qual[GT_Qual.fixed != 1]
            
            fileName= args.o + '_' + p + '.csv'
            # if file does not exist write header 
            if not os.path.isfile(fileName):
                GT_Qual.to_csv(fileName, mode = 'w', index=False)
            else: # else it exists so append without writing the header
                GT_Qual.to_csv(fileName, mode = 'a',header=False, index=False)
        print('[DONE]')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Format Data for Regression Analysis')
    parser.add_argument('-i', help = 'input Directory')
    parser.add_argument('-o', help= 'output directory')
    args = parser.parse_args()
    main(args)
    
    
    #for f in `seq 1 21`; do echo $f; java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_Chr${f}.4bed_filtered.vcf.gz  CHROM POS "GEN[*].GT"   | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d'| tail -n +2 > /Users/luke/Documents/Regression/CHR${f}.Genotypes.txt; done