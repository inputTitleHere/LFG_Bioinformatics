from numpy.lib.npyio import save
import pandas as pd
import numpy as np
import os
from statsmodels.stats import multitest as mt
from pathlib import Path
#FDR coding?

def location_list(gene_name):
    root=Path(__file__).parents[1]
    #path=os.path.join(root, 'Test_Output',gene_name)
    path=os.path.join(root, 'Output_Results')
    files_list=os.listdir(path)
    files_path_list=[os.path.join(path, loc) for loc in files_list if gene_name.upper() in loc.upper()]
    return files_path_list

def open_dataframe(location):
    df=pd.read_csv(location,sep=',')
    return df

def absolutify(df):
    mean_diff=df['mean_diff'].copy()
    count=0
    pos_neg=[]
    for i in mean_diff:
        count+=1
        if i>0:
            pos_neg.append('+')
        elif i==0:pos_neg.append('0')
        elif i<0:pos_neg.append('-')
    pos_neg_s=pd.Series(pos_neg, name='pos_neg')
    df=pd.concat([df,pos_neg_s],axis=1)
    df['mean_diff']=df['mean_diff'].map(lambda x:abs(x))
    df=df.sort_values(by=['mean_diff'], ascending=False)
    return df

def save_to_csv_2(df,gene_name,name='Output_Absolutify', index=False):
    loc=Path(__file__).parents[1]
    #Add '_abs' to namespace befor extension information.
    loc=os.path.join(loc, 'Output_Absolutify',gene_name, f'{name}_abs.csv')
    df.to_csv(loc, sep=',',index=index)
    return


if __name__=='__main__':
    script_path=Path(__file__).parents[1]
    #gene_name_list=os.listdir(os.path.join(script_path, 'Test_Output'))
    wnt_path='C://Users//c//Desktop//LFG_Bioinformatics//Target_list//wnt_target_genes_2.txt'
    wnt_target_genes=pd.read_csv(wnt_path, sep='\t', encoding='utf-8', names=['Gene_symbol','Ensembl_ID'])
    gene_name_list=list(pd.Series(wnt_target_genes['Gene_symbol']))
    for gene_name in gene_name_list:
        # If folder with gene name does not exist, create folder within Test_Output_2
        if not os.path.exists(os.path.join(Path(__file__).parents[1],'Output_Absolutify',gene_name)):
            Path(os.path.join(script_path,'Output_Absolutify',gene_name)).mkdir(parents=True, exist_ok=True)
        file_list=location_list(gene_name)
        print(f'Processing {gene_name}')
        for file in file_list:
            cancer_name=file[file.rfind("\\")+1:file.rfind('_')]
            df=open_dataframe(file)
            df=absolutify(df)
            save_to_csv_2(df, gene_name=gene_name,name=f'{cancer_name}_{gene_name}')

