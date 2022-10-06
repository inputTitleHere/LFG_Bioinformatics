from typing import ClassVar
import mygene as myg
from numpy.core.fromnumeric import sort
import pandas as pd
import os
'''
target_gene={
    'CCN2':['ENSG00000118523'], 'CCN1':['ENSG00000142871'],'CCND1':['ENSG00000110092'],'MAPK1':['ENSG00000100030'],'HES1':['ENSG00000114315'],'NRARP':['ENSG00000198435'],'HK1':['ENSG00000156515'],'GPI':['ENSG00000105220'],
    'NKD1':['ENSG00000140807'],'Axin2':['ENSG00000168646'],'APC':['ENSG00000134982']}

#print(list(target_gene.keys()))
#print(list(target_gene.values()))
names_lst=[i[0] for i in target_gene.values()]
'''

def open_fpkm(cancer,gene_list,dic,sort_symbols):
    pos='C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/FPKM_Data'
    files=os.listdir(pos)
    locs=[i for i in files if cancer in i]
    path=os.path.join(pos,locs[0])
    df=pd.read_csv(path, sep='\t',index_col=False)
    df['Ensembl_ID']=df['Ensembl_ID'].map(lambda x:x[:x.find('.')]).copy()
    df2=df.loc[df['Ensembl_ID'].isin(gene_list)].copy()
    df2['Ensembl_ID']=df2['Ensembl_ID'].map(lambda x:dic[x]).copy()
    df2=df2.set_index('Ensembl_ID')
    df2=df2.loc[sort_symbols].copy()
    df2.to_csv(f'C:/Users/c/Desktop/t/{cancer}.csv')



df=pd.read_csv('C:/Users/c/Desktop/t/div_15.csv', sep=',',index_col=False)
gene_sym=list(df['Gene_symbol'])
print(gene_sym)

convert={}



mg=myg.MyGeneInfo()
rslt=mg.querymany(gene_sym,scopes='symbol',fields='ensembl.gene',species='human')
for items in rslt:
    try:
        convert[items['ensembl']['gene']]=items['query']
    except TypeError:
        convert[items['ensembl'][1]['gene']]=items['query']
print(convert)
gene_list=list(convert.keys())
print(gene_list)
cancer_list=['UCEC','BRCA']

for c in cancer_list:
    open_fpkm(c,gene_list=gene_list,dic=convert,sort_symbols=gene_sym)