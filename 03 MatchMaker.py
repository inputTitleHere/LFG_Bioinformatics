import os
from os.path import basename
import pandas as pd
import numpy as np
from pathlib import Path

def open_df(path):
    '''Open dataframe with index inclueded as default. Index will be used to determine the rank of each genes.'''
    df=pd.read_csv(path, sep=',')
    return df

def save_to_csv(path,df,name,index=False):
    '''Save dataframe to csv. Place it at Output_3 folder.'''
    name_=f'{name[:name.rfind(".")]}_meancut.csv'
    loc=os.path.join(path, name_)
    df.to_csv(loc, sep=',',index=index)

if __name__=='__main__':
    path=Path(__file__).parents[1]
    #Variable(Gene name) goes here
    #gene_list=os.listdir(os.path.join(path, 'Test_Output_2'))
    gene_list=os.listdir(os.path.join(path, 'Output_Absolutify'))
    for ref_gene_name in gene_list:
        print(ref_gene_name)
        path_read=os.path.join(path,'Output_Absolutify',ref_gene_name)
        path_save=os.path.join(path, 'Output_MatchMaker')
        if not os.path.exists(path_save):
            Path(path_save).mkdir(parents=True, exist_ok=True)
        file_list=os.listdir(path_read)
        print(file_list)
        #Cancer name list used for getting columns.
        cancer_list=[]
        #Create dataframe to save all combined information of a single pathway target gene.
        gene_df=pd.DataFrame(columns=['gene_info'])
        for file in file_list:
            file_path=os.path.join(path_read, file)
            #basename() gets file name from a string of path.
            base_name=os.path.basename(file_path)[:file_path.rfind('.')]
            #Get cancer name from file name. 
            cancer_name=base_name[:base_name.find("_")]
            cancer_list.append(cancer_name)
            print(f'==========================={cancer_name}============================')
            df=open_df(path=file_path)
            #df=df.loc[df['mean_diff']>0.1]
            for i, s in df.iterrows():
                #Merged gene names into a single string for convenience
                gene_name=f'{s["symbol"]}_[{s["splice_type"]}]_"{s["exons"]}"({s["from_exon"]}->{s["to_exon"]})'
                #Merged PSI values into a single string. (Split to different cells.0310)
                psi_info=f"{s['mean_diff']:.6f}({s['pos_neg']})"
                if gene_name not in gene_df['gene_info'].values:
                    gene_df=gene_df.append({'gene_info':gene_name, f'{cancer_name}':psi_info}, ignore_index=True)
                elif gene_name in gene_df.values:
                    gene_df.loc[gene_df['gene_info']==gene_name,f'{cancer_name}']=psi_info
        print('__Calculating Dataframe___')
        #Code for getting mean values of each PSI.
        #First check if all values in rows are positive. 
        gene_df[cancer_list]=gene_df[cancer_list].applymap(lambda x: float(f"{x[x.find('(')+1:x.rfind(')')]}{x[:x.find('(')]}"), na_action='ignore')

        #Count cancer samples with value
        gene_df['cancer_count']=len(cancer_list)-gene_df.isna().sum(axis=1)
        #Calculate delta mean
        delta_mean=[]
        for i, s in gene_df.iterrows():
            #cancer_list dataframe above is saved as float. dtype below must be specified as float.
            #Get delta mean value excluding nan data types
            delta_mean_value=np.nanmean(s[cancer_list],dtype='float')
            delta_mean.append(delta_mean_value)
        delta_mean_sereis=pd.Series(delta_mean)
        gene_df['delta_mean']=delta_mean_sereis
        gene_df['eval_score']=abs(gene_df['delta_mean']*gene_df['cancer_count'])
        #Sort dataframe by cancer count
        gene_df.sort_values(by='eval_score', ascending=False,inplace=True)
        #gene_df=gene_df.loc[gene_df['cancer_count']>=len(cancer_list)]
        gene_df=gene_df.loc[abs(gene_df['delta_mean'])>0.1]
        print(f'Saving file of {ref_gene_name}\n')
        save_to_csv(path_save, gene_df, name=f'{ref_gene_name}_merged',index=False)
 