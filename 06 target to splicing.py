from operator import is_not, truediv
from os import path
from numpy.core.fromnumeric import compress
from numpy.lib.arraysetops import intersect1d
import pandas as pd
import os
from pathlib import Path
# This script is for getting pathway-target-gene (vs) spliced gene(with all splice variants)

def open_FPKM(path, ensembl_id, gene_symbol):
    '''Return dataframe containing gene of interest. Input ensembl ID as ensembl_id'''
    cancer_name=path[path.find('-')+1:path.find('.ht')]
    ensembl_list=[*ensembl_id]
    print(f"Processing _______[__{cancer_name}__]________")
    fpkm_df=pd.read_csv(path, sep='\t',compression='gzip', encoding='utf-8')
    fpkm_df['Ensembl_ID']=fpkm_df['Ensembl_ID'].map(lambda x:x[:x.find('.')]).copy()
    interest_df=fpkm_df.loc[fpkm_df['Ensembl_ID'].isin(ensembl_list)]
    name_row=list(interest_df.columns)
    #Ensembl_ID.
    title_column=[name_row[0]]
    canc_row=[]
    #norm_row=[]
    #Patient ID
    for names in name_row:
        if '-11A' in names or '-11B' in names:
            #norm_row.append(names)
            continue
        if '-01A' in names or '-01B' in names or '-01C' in names or '-02A' in names:
            canc_row.append(names)
            continue
    #Asterisk * is used to unzip lists.
    sorted_row=[*title_column,*canc_row]
    interest_df=interest_df[sorted_row]
    column_names=interest_df.columns
    column_names=map(lambda x:x[:x.rfind('-')] if '-' in x else x, column_names)
    column_names=map(lambda x:x.replace('-','_') if '-' in x else x, column_names)
    interest_df.columns=column_names
    interest_df.set_index(keys='Ensembl_ID', drop=True, inplace=True)
    interest_df=interest_df.rename_axis(None)
    interest_df=interest_df.transpose()
    interest_df.sort_values(by=ensembl_list, ascending=False, inplace=True, ignore_index=False)
    interest_df.columns=[gene_symbol]
    return interest_df, cancer_name, list(interest_df.index)

def open_psi(cancer_name, gene_symbol='MARK3'):
    location=os.path.join(Path(__file__).parents[1], 'Test_data','PSI_Data')
    #one line code to find splice data contaning correct cancer names. Hope it works
    psi_file_list=[t for t in os.listdir(location) if cancer_name in t]
    psi_file=psi_file_list[0]
    file_path=os.path.join(location, psi_file)
    psi_df=pd.read_csv(file_path, compression='zip', sep='\t',encoding='utf-8')
    psi_df=psi_df.loc[psi_df['symbol']==gene_symbol]
    psi_columns=list(psi_df.columns)
    canc=[]
    # Get the list of cancer sample names.
    for name in psi_columns:
        if 'TCGA' in name and '_Norm' in name:
            #Omitting Normal list.
            continue
        if 'TCGA' in name and not '_Norm' in name:
            canc.append(name)
    # Only collect names list and cancer samples. Omit normal samples. 
    psi_info=['symbol','splice_type','exons','from_exon','to_exon']
    name_fin_lst=[*psi_info, *canc]
    psi_df=psi_df[name_fin_lst].copy()
    gene_info=[]
    for _,s in psi_df.iterrows():
        gene_info.append(f"{s['symbol']}_[{s['splice_type']}]_\"{s['exons']}\"({s['from_exon']}->{s['to_exon']})")
    psi_df=psi_df[canc].copy()
    psi_df.insert(loc=0, column='gene_info', value=gene_info)
    psi_df.set_index(keys='gene_info', drop=True, inplace=True)
    psi_df=psi_df.rename_axis(None)
    psi_df=psi_df.transpose()
    return psi_df, list(psi_df.index)


def match_patients(fpkm_patient_list, psi_patient_list):
    matching_list=[i for i in fpkm_patient_list if i in psi_patient_list]
    match_df=pd.DataFrame(matching_list, columns=['ids'])
    #Going through unnessacerily long codes to completely remove duplicate samepls
    tdf=match_df.duplicated()
    tdf2=tdf.loc[tdf==True]
    dup=list(tdf2.index) # Contains numeric coordinates of indexes
    dup_df=match_df.iloc[dup]
    duplist=list(dup_df['ids'])
    matching_both=match_df.loc[~match_df['ids'].isin(duplist)]
    print(f"FPKM patient sample count = {len(fpkm_patient_list)}\nPSI patient sample count = {len(psi_patient_list)}\nMathicng list = {len(matching_both)}")
    return list(matching_both['ids'])

def select_percentage(fpkm_psi_df, cutoff_percentage=0.2, exons=None, from_e=None, to_e=None):
    '''Input cutoff threshold of either 1~100 or 0~1. 1 will be counted as 1%'''
    #Check if percentage or float value(100, 0.1 ect)
    if cutoff_percentage<100 and cutoff_percentage>=1:
        cutoff_percentage=cutoff_percentage/100
    else:
        pass
    #Check if exons are inputted
    if exons is not None and from_e is not None and to_e is not None:
        from_e=float(from_e)
        to_e=float(to_e)
        exon_interest=f'"{exons}"({from_e}->{to_e})'
        splice_interest=[c for c in fpkm_psi_df.columns if exon_interest in c]
        fpkm_psi_df=fpkm_psi_df[splice_interest].copy()
    else:
        pass
    index_count=len(fpkm_psi_df.index)
    num_index=round(index_count*cutoff_percentage)
    topside_df=fpkm_psi_df.iloc[:num_index]
    bottomside_df=fpkm_psi_df.iloc[index_count-num_index:]
    print(topside_df)
    print(bottomside_df)
    return topside_df, bottomside_df

def save_to_csv(df, gene_symbol,name='(Cancer_name)_(gene_name)_(gene_of_interest)', index=False):
    loc=Path(__file__).parents[1]
    path_save=os.path.join(loc, 'Output_FPKM_PSI_match',gene_symbol)
    if not os.path.exists(path_save):
        Path(path_save).mkdir(parents=True, exist_ok=True)
    loc=os.path.join(path_save, f'{name}.csv')
    df.to_csv(loc, sep=',',index=index)
    return


if __name__=='__main__':
    fpkm_target=[]
    psi_target=[]

    fpkm_path=os.path.join(Path(__file__).parents[1], 'Test_data','FPKM_Data')
    file_list=os.listdir(fpkm_path)
    #=============#
    target_gene_2={'fra-1':['ENSG00000175592'], 'PPARD':['ENSG00000112033'], 'PROM1':['ENSG00000007062']}
    #target_gene_2={'EGFR':['ENSG00000146648'],'14-3-3_Sigma':['ENSG00000175793']}
    gene_of_interest=['CD47','IRF3','RALGAPA1','SYTL2']
    #=============#
    file_name=[]
    cancer_name_list=[]
    for file in file_list:
        if '.gz' in file:
            file_name.append(file)

    # open RNAseq(FPKM)
    for gene_symbol, ensembl_id in target_gene_2.items():
        for fpkm_name in file_name:
            for gene_of_interest_single in gene_of_interest:
                print(f'========[{gene_symbol}]========')
                loc=os.path.join(fpkm_path, fpkm_name)
                # CCN1 ='ENSG00000142871'
                fpkm_df, cancer_name, fpkm_patient_ids=open_FPKM(loc, ensembl_id=ensembl_id, gene_symbol=gene_symbol)
                psi_df, psi_patient_ids=open_psi(cancer_name=cancer_name, gene_symbol=gene_of_interest_single)
                print(f'{cancer_name}, {gene_symbol}, {gene_of_interest_single}')
                #print(fpkm_df)
                #print(psi_df)
                match_list=match_patients(fpkm_patient_list=fpkm_patient_ids, psi_patient_list=psi_patient_ids)
                fpkm_psi_df=pd.concat([fpkm_df.loc[match_list],psi_df.loc[match_list]],axis=1) 
                namespace=f'{cancer_name}_{gene_symbol}_{gene_of_interest_single}'
                fpkm_psi_df.index.set_names('patient_id',inplace=True)
                #topside_df, bottomside_df=select_percentage(fpkm_psi_df=fpkm_psi_df,cutoff_percentage=20, exons=17, from_e=16, to_e=18)
                save_to_csv(fpkm_psi_df, gene_symbol=gene_symbol, name=namespace,index=True)
                print(f"Finished Processing {cancer_name}\n\n")

print("All process finished")
