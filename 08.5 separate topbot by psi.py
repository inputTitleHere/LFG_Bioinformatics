import os
from numpy.lib.shape_base import expand_dims
import pandas as pd
from pathlib import Path
from scipy.stats import ttest_ind as tt

cancer_name_list=os.listdir('C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/PSI_Data')
cancer_name_list=list(map(lambda x:x[x.rfind('_')+1:x.find('.')],cancer_name_list))

fpkm_psi_merge_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Output_FPKM_PSI_match'

fpkm_target_gene_list=['fra-1','PPARD','PROM1']
symbol_exon_info_df=pd.read_csv('C:/Users/c/Desktop/LFG_Bioinformatics/Target_list/wnt_symbol_exon_info.txt', sep='\t')

#===#
total_ccn1_list=[]
total_ccn2_list=[]
#===#

for idx,rows in symbol_exon_info_df.iterrows():
    list_row=list(rows)
    splice_gene=list_row[0]
    exon=list_row[1]
    from_exon=float(list_row[2])
    to_exon=float(list_row[3])
    splice_info=f'"{exon}"({from_exon}->{to_exon})'

for names in cancer_name_list:
    print(f'___Processing {names}___')
    file_path_ccn1=f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN1/{names}_CCN1.csv'
    file_path_ccn2=f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN2/{names}_CCN2.csv'
    # Read csv from file
    df_ccn1=pd.read_csv(file_path_ccn1, sep=',', encoding='utf-8')
    df_ccn2=pd.read_csv(file_path_ccn2, sep=',', encoding='utf-8')
    # Find if data contains certain splice type
    splice_type=[i for i in df_ccn1.columns if splice_info in i]
    if not splice_type:
        # If certain splice type not exists, skip to the next cancer
        continue
    # Get dataframe where splice value is not null(data exists)
    df_ccn1=df_ccn1[~df_ccn1[f'{splice_type[0]}'].isnull()]
    df_ccn2=df_ccn2[~df_ccn2[f'{splice_type[0]}'].isnull()]
    # splice by splice type 
    df_ccn1.sort_values(splice_type[0], inplace=True, ascending=False)
    df_ccn2.sort_values(splice_type[0], inplace=True, ascending=False)
    df_ccn1.set_index('patient_id', drop=True, inplace=True)
    df_ccn2.set_index('patient_id', drop=True, inplace=True)
    indexes=list(df_ccn1.index)
    df_ccn2=df_ccn2.loc[indexes]
    #===#
    # Getting 20% of each end
    locus=round(len(df_ccn1.index)*0.2)
    #===#
    #=========================================================================#
    ccn1_df=df_ccn1['CCN1']
    ccn1_df=ccn1_df.to_frame()
    #ccn1_df.rename({'CCN1':f'{names}_CCN1'}, inplace=True)
    # top side is where mark3 psi is near 1 (high)
    top_ccn1_df=ccn1_df.iloc[:locus].copy()
    top_ccn1_df.reset_index(drop=True, inplace=True)
    top_ccn1_df.rename(columns={'CCN1':f'top_{names}_CCN1'}, inplace=True)
    # bot side is where mark3 psi is low (near 0)
    bot_ccn1_df=ccn1_df.iloc[len(df_ccn1.index)-locus:].copy()
    bot_ccn1_df.reset_index(drop=True, inplace=True)
    bot_ccn1_df.rename(columns={'CCN1':f'bot_{names}_CCN1'}, inplace=True)
    #=# appending to total list for concat
    total_ccn1_list.append(top_ccn1_df)
    total_ccn1_list.append(bot_ccn1_df)
    tstat, pval=tt(top_ccn1_df, bot_ccn1_df)
    print(f'CCN1 pvalue = {pval}') 
    #==========================================================================#
    ccn2_df=df_ccn2['CCN2']
    ccn2_df=ccn2_df.to_frame()
    #ccn2_df.rename({'CCN2':f'{names}_CCN2'}, inplace=True)
    top_ccn2_df=ccn2_df.iloc[:locus].copy()
    top_ccn2_df.reset_index(drop=True, inplace=True)
    top_ccn2_df.rename(columns={'CCN2':f'top_{names}_CCN2'}, inplace=True)
    # bot side is where mark3 psi is low (near 0)
    bot_ccn2_df=ccn2_df.iloc[len(df_ccn2.index)-locus:].copy()
    bot_ccn2_df.reset_index(drop=True, inplace=True)
    bot_ccn2_df.rename(columns={'CCN2':f'bot_{names}_CCN2'}, inplace=True)
    #=# appending to total list for concat
    total_ccn2_list.append(top_ccn2_df)
    total_ccn2_list.append(bot_ccn2_df)
    tstat, pval=tt(top_ccn2_df, bot_ccn2_df)
    print(f'CCN2 pvalue = {pval}') 
    #==========================================================================#
    splice_columns=df_ccn1.columns
    splice_of_interest=[i for i in splice_columns if splice_info in i]
    splice_df=df_ccn1[splice_of_interest[0]]
    splice_df=splice_df.to_frame()
    fin_df=pd.concat([splice_df, ccn1_df, ccn2_df], axis=1)
    
    #fin_df.to_csv(f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_6/{names}_MARK3_{exons}({from_e}_{to_e}).csv',index=False, encoding='utf-8')
    print('\n\n')

ccn1_total=pd.concat(total_ccn1_list, axis=1)
#ccn1_total.to_csv(f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_6/MARK3_{exons}({from_e}_{to_e}_CCN1).csv',index=False, encoding='utf-8')
ccn2_total=pd.concat(total_ccn2_list, axis=1)
#ccn2_total.to_csv(f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_6/MARK3_{exons}({from_e}_{to_e}_CCN2).csv',index=False, encoding='utf-8')