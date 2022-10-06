from os import path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from statistics import mean
from pathlib import Path
from scipy import stats


def get_FPKM_loc_list():
    '''Open FPKM data folder and return a list of each data path.'''
    root=Path(__file__).parents[1]
    #root='LFG_BIOINFORMATICS'
    #Find folder 'Test data' and 'FPKM data' folder inside. 
    path=os.path.join(root, 'Test_data','FPKM_Data')
    files=os.listdir(path)
    #fpkm_loc=[loc for loc in files if 'fpkm' in loc.lower()]
    fpkm_path_list=[os.path.join(path, loc) for loc in files if loc.endswith('.gz')]
    return fpkm_path_list

def open_fpkm(fpkm_loc):
    '''Open RNA Expression Data(FPKM) into a pandas dataframe'''
    fpkm_path=fpkm_loc
    # Open Dataframe here.
    df=pd.read_csv(fpkm_path, sep='\t', encoding='utf-8')
    # /Open Dataframe here

    new_ensembl_frame=df.loc[:,'Ensembl_ID'].map(lambda x : x[:x.index('.')])
    df['Ensembl_ID']=new_ensembl_frame
    sorted_row=[]
    name_row=list(df.columns)
    #Ensembl_ID.
    title_column=[name_row[0]]
    norm_row=[]
    canc_row=[]
    #Patient ID
    for names in name_row:
        if '-11A' in names or '-11B' in names:
            norm_row.append(names)
            continue
        if '-01A' in names or '-01B' in names or '-01C' in names or '-02A' in names:
            canc_row.append(names)
            continue
    #Asterisk * is used to unzip lists.
    sorted_row=[*title_column,*norm_row,*canc_row]
    df=df[sorted_row]
    df.set_index('Ensembl_ID', inplace=True)
    return df, len(title_column), len(norm_row), len(canc_row)

def open_psi(cancer_name):
    root=Path(__file__).parents[1]
    path=os.path.join(root, 'Test_data','PSI_Data')
    files=os.listdir(path)
    for loc in files:
        # Open Dataframe with matching Cancer name. 
        if loc.endswith('.txt') and 'psi' in loc.lower() and cancer_name.upper() in loc.upper():
            psi_path=os.path.join(path, loc)
            df=pd.read_csv(psi_path, sep='\t', encoding='utf-8')
        elif loc.endswith('.zip') and 'psi' in loc.lower() and cancer_name.upper() in loc.upper():
            psi_path=os.path.join(path, loc)
            df=pd.read_csv(psi_path, compression='zip',sep='\t', encoding='utf-8')
    # Fixed list of names list for all data downloaded from TCGA SpliceSeq.
    names_list=['symbol','splice_type','exons','from_exon','to_exon']
    name_row=list(df.columns)
    canc=[]
    # Get the list of cancer sample names.
    for name in name_row:
        if 'TCGA' in name and '_Norm' in name:
            #Omitting Normal list.
            continue
        if 'TCGA' in name and not '_Norm' in name:
            canc.append(name)
    # Only collect names list and cancer samples. Omit normal samples. 
    name_fin_lst=[*names_list, *canc]
    df=df[name_fin_lst].copy()
    return df

def exctract_matching(main_df, gene_dic):
    '''Get dataframe of a certain gene.'''
    # gene_dic is a dictionary with values as list.(incase multiple ensembl IDs for single gene symbol. in build v01.03/12/21 the dictionary is fixed.)
    names_lst=[i[0] for i in gene_dic.values()]
    # chage everytime // current gene =[CCND1]
    # get rows that contain genes of interest into a dataframe
    df=main_df[main_df.index.isin(names_lst)].copy()
    return df

def get_patient_ID(row_df, range_cutoff=0.20):
    '''Default Percentage is 20%'''
    indexing_flag=False
    percentage_flag=True
    #Change Values Here!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #indexing flag is not to be used in a full run. Only for probing sample percentage manually.
    if indexing_flag:
        low_index=list(row_df[row_df<16].keys())
        high_index=list(row_df[row_df>28.7].keys())
        return [low_index, high_index]
    # Using indexes as i am getting lazy
    elif percentage_flag:
        index_list=list(row_df.keys())
        low_level=index_list[:round(len(index_list)*range_cutoff)]
        high_level=index_list[len(index_list)-round(len(index_list)*range_cutoff):]
        return [low_level, high_level]

def match_psi(psi_df, matching_list):
    '''Match PSI and patient IDs.'''
    tot_list1=matching_list[0]
    tot_list2=matching_list[1]
    #map items within each list. manipulating patient IDs
    tot_list1=list(map(lambda x:x[:x.rfind('-')], tot_list1))
    tot_list1=list(map(lambda x:x.replace('-','_'),tot_list1))
    tot_list2=list(map(lambda x:x[:x.rfind('-')], tot_list2))
    tot_list2=list(map(lambda x:x.replace('-','_'),tot_list2))
    '''
    t_list=[*matching_list[0],*matching_list[1]]
    t_list=list(map(lambda x:x[:x.rfind('-')], t_list))
    t_list=list(map(lambda x:x.replace('-','_'),t_list))
    '''
    psi_list=list(psi_df.columns)
    #===================================================#
    #                                                   #
    #===================================================#
    #name_list=['symbol','splice_type','exons','from_exon','to_exon',*t_list]
    #!!!! When deleting stuff from lists, make a copy of it. !!!!
    for items in tot_list1[:]:
        if items not in psi_list:
            tot_list1.remove(items)
    for items in tot_list2[:]:
        if items not in psi_list:
            tot_list2.remove(items)
    print(f'Length of lower list = {len(tot_list1)}, Length of higher list = {len(tot_list2)}.\n')
    basic_column_names=['symbol','splice_type','exons','from_exon','to_exon']
    fin_name_list=[*basic_column_names,*tot_list1,*tot_list2]
    df=psi_df[fin_name_list]
    return df, [basic_column_names,tot_list1,tot_list2]

def calculate_pval(df, column_info):
    count=0
    pval_res=[]
    mean_low_list=[]
    mean_high_list=[]
    diff_list=[]
    lower_count={'Lower count':[]}
    Higher_count={'Higher count':[]}
    for index, series in df.iterrows():
        count+=1
        lower_not_NaN_count=0
        higher_not_NAN_count=0
        low=series[column_info[1]].dropna()
        high=series[column_info[2]].dropna()
        try:
            result=stats.ttest_ind(low, high, equal_var=False)
            pval_res.append(result.pvalue)
        except ZeroDivisionError:
            #print(f'Count={count}, Mean={mean_val}')
            pval_res.append(1)
        #print(result.pvalue)
        mean_low=mean(low)
        mean_high=mean(high)
        mean_low_list.append(mean_low)
        mean_high_list.append(mean_high)
        diff=mean_high-mean_low
        diff_list.append(diff)
    p_value_s=pd.Series(pval_res, name='p_value')
    mean_low_s=pd.Series(mean_low_list,name='mean_low')
    mean_high_s=pd.Series(mean_high_list,name='mean_high')
    mean_diff_s=pd.Series(diff_list,name='mean_diff')
    n_df=pd.concat([df, p_value_s, mean_low_s, mean_high_s, mean_diff_s], axis=1)
    return n_df

def remove_patients(calculated_df):
    '''Remove patients ID from dataframe and sort by Delta'''
    names=['symbol','splice_type','exons','from_exon','to_exon','p_value','mean_low','mean_high','mean_diff']
    calculated_df=calculated_df[names].copy()
    calculated_df=calculated_df[calculated_df['p_value']<0.05]
    calculated_df=calculated_df.replace(':',';',regex=True)
    cal_df=calculated_df.sort_values(by=['mean_diff'], ascending=False)
    return cal_df

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

def save_to_csv(df, gene_name='(Gene_name)',name='(Cancer_name)_(gene_name)', index=False):
    loc=Path(__file__).parents[1]
    path_save=os.path.join(loc, 'Test_Output', gene_name)
    if not os.path.exists(path_save):
        Path(path_save).mkdir(parents=True, exist_ok=True)
    loc=os.path.join(loc, 'Test_Output', gene_name, f'{name}.csv')
    df.to_csv(loc, sep=',',index=index)
    return

#==========================Target Genes=============================#
'''
target_gene={'CCN2':['ENSG00000118523'], 'CCN1':['ENSG00000142871'],'CCND1':['ENSG00000110092'],'HES1':['ENSG00000114315'],'NRARP':['ENSG00000198435'],'HK1':['ENSG00000156515'],'GPI':['ENSG00000105220'],'NKD1':['ENSG00000140807'],'Axin2':['ENSG00000168646'],'APC':['ENSG00000134982'],'RBFOX2':['ENSG00000100320']}
'''
#'GSK3b':['ENSG00000082701'],'DVL1':['ENSG00000107404'],'LEF1':['ENSG00000138795']
target_gene_2={'CCN2':['ENSG00000118523'], 'CCN1':['ENSG00000142871']}
ensembl_loc=os.path.join(Path(__file__).parents[1], 'Target_to_Ensembl')
keyword='hippo'
target_gene={}
flagish=False
if flagish:
    for file in os.listdir(ensembl_loc):
        if keyword.lower() in file.lower():
            with open(os.path.join(ensembl_loc, file), mode='r', encoding='utf-8') as f:
                for ensembl_ids in f.readlines():
                    target_gene[ensembl_ids[:ensembl_ids.find('\t')]]=[ensembl_ids[ensembl_ids.find('\t')+1:ensembl_ids.find('\n')]]

#Set using gene list here.
gene_space=target_gene_2

#Set percentage of top and bottom here.
range_cutoff=0.15
#==========================MAIN SPACE=============================#

if __name__=='__main__':
    f_p=True
    #Get list of cancer fpkm data.
    # Input data for FPKM is retrieved from https://xenabrowser.net/datapages/ GDC TCGA{Cancer name}/gene expression RNAseq/HTSeq-FPKM(not FPKM-UQ)
    fpkm_file_loc_list=get_FPKM_loc_list()
    for fpkm_loc in fpkm_file_loc_list:
        #Get cancer name as string.
        cancer_name=fpkm_loc[fpkm_loc.find('-')+1:fpkm_loc.find('.')]
        print(f'____________Processing {cancer_name}____________')
        # With cancer name retrieved from above, open the matching PSI dataframe.
        # PSI data is retrived from https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIdownload.jsp by selecting cancer type, ES only(for now) and percentage of samples with PSI value of 75.
        psi_df=open_psi(cancer_name=cancer_name)
        #base name is file name.
        basename=os.path.basename(fpkm_loc)
        #FPKM dataframe for single cancer, and column info as len_list
        fpkm_df, *len_list=open_fpkm(fpkm_loc)
        fpkm_n_df=exctract_matching(fpkm_df, gene_space)
        print(fpkm_n_df) #row count = number of target genes. 
        for index, series in fpkm_n_df.iterrows():
            print(f'Ensembl ID = {index}')
            # index contains Ensmbl IDs of interest.
            # k=Key, v=Value of dictionary {key:value} -> {'gene symbol':['ensembl_ID']}
            gene_key=[k for k,v in gene_space.items() if v[0]==index]
            #return string from list. The list will contain one item only anyways.
            gene_name=gene_key[0]
            print(f'Gene Symbol = {gene_name}')
            #for series(row)in fpkm dataframe: get name and cancer columns only(remove normal)
            # Remove name column and normal column from series by indexing iloc
            row=series.iloc[(len_list[0]+len_list[1]):].copy()
            row=row.sort_values()
            # Get list of patients for top and bottom(20% for now)
            id_list=get_patient_ID(row,range_cutoff=range_cutoff)
            # extract df of top/bottom % into matched_df. column_info includes basic column names and top/bottom patient IDs. 
            matched_df, column_info=match_psi(psi_df, id_list)
            # Calculate pvalue of PSI for matching genes. 
            calculated_df=calculate_pval(matched_df, column_info)
            #remove patient IDs and extract information only. 
            extracted_df=remove_patients(calculated_df)
            #Save to csv by the name (Cancer_name)_(Gene_name)
            save_to_csv(extracted_df, gene_name=gene_name, name=f'{cancer_name}_{gene_name}',index=False)