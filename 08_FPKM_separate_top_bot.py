from operator import index
from sys import path
import pandas as pd
from pathlib import Path
import os
from pandas.core.reshape.merge import merge
from scipy.stats import ttest_ind as tt


# defs goes here

class Separate_top_bottom:

  def __init__(self):
    self.list_of_dfs=[]
    self.list_of_heatmapping=[]
    self.top_df=None
    self.bot_df=None
    self.tt_res={}
    self.final_df=None
    self.current_fpkm=''
    self.current_splicing=''
    self.current_cancer=''
    print('===== Initializing Class =====')

  def Separation_boxploting(self, location, fpkm_gene, splice_of_interest):
    self.current_fpkm=fpkm_gene
    self.current_splicing=splice_of_interest
    df=pd.read_csv(location, sep=',',encoding='utf-8')
    splice_type=[i for i in df.columns if splice_of_interest in i]
    if not splice_type:
      # If certain splice type not exists, skip to the next cancer
      return
    df=df[~df[splice_of_interest].isnull()]
    df.sort_values(splice_of_interest,inplace=True,ascending=False)
    df.set_index('patient_id',drop=True, inplace=True)
    list_of_index=list(df.index)
    locus=round(len(list_of_index)*0.2)
    series_fpkm=df[fpkm_gene].to_frame()
    # top df starts here
    self.top_df=series_fpkm.iloc[:locus].copy()
    self.top_df.reset_index(drop=True, inplace=True)
    self.top_df.rename(columns={fpkm_gene:f'T_{cancer_name}_{fpkm_gene}'},inplace=True)
    # bot df starts here
    self.bot_df=series_fpkm.iloc[len(df.index)-locus:].copy()
    self.bot_df.reset_index(drop=True, inplace=True)
    self.bot_df.rename(columns={fpkm_gene:f'B_{cancer_name}_{fpkm_gene}'},inplace=True)
    # T-test for stars :)
    tstat, pval=tt(self.top_df, self.bot_df)
    self.tt_res[f'{cancer_name}_{fpkm_gene}_{splice_of_interest}']=pval[0]
    # Throw top bot df to the list for future merging
    self.list_of_dfs.extend((self.top_df, self.bot_df))
    return

  def Separate_heatmapping(self, location_merge_loc, fpkm_gene_list, splice_gene, exon, splice_of_interest,cancer_type):
    self.current_splicing=splice_gene
    root_loc_list=os.listdir(location_merge_loc)
    storage_dic={}
    folder_target=[os.path.join(location_merge_loc, i) for i in root_loc_list if i in fpkm_gene_list]
    # First index on fpkm_gene_list will be used as reference for other dataframes.
    for location in folder_target:
      folder_list=os.listdir(location)
      for fpkm_gene in fpkm_gene_list:
        name=f'{cancer_type}_{fpkm_gene}_{splice_gene}.csv'
        if name in folder_list:
          df=pd.read_csv(os.path.join(location, name))
          storage_dic[f'{fpkm_gene}']=df.set_index('patient_id',drop=True)
    index_list=list(storage_dic[fpkm_gene_list[0]].index)
    n=0
    merge_list=[]
    for k,v in storage_dic.items():
      storage_dic[k]=v.loc[index_list]

      if n==0:
        merge_list.append(v.get(splice_of_interest))
        merge_list.append(v[fpkm_gene_list[n]].to_frame())
      else:
        merge_list.append(v[fpkm_gene_list[n]].to_frame())
      n+=1
    merge_df=pd.concat(merge_list, axis=1)
    self.final_df=merge_df
    return

  def Save_file_boxplot(self, path):
    save_df=pd.concat(self.list_of_dfs, axis=1)
    s=self.current_splicing.replace(':',';')
    splice_name=s[:s.find('_')+1]+s[s.find('\"')+1:s.rfind('\"')]
    save_path=os.path.join(path, 'Boxplotting',f'TopBot_{self.current_fpkm}_{splice_name}.csv')
    save_df.to_csv(path_or_buf=save_path, index=False)

    for k,v in self.tt_res.items():
      v=str(v)
      v=v.replace(' ','')
      if 'e' in v:
        self.tt_res[k]=v[:v.find('e')]+' '+v[v.find('e'):]
    pvals_df=pd.DataFrame.from_dict(self.tt_res, orient='index')
    pvals_df.to_csv(path_or_buf=os.path.join(path, 'Boxplotting','pvalue_results.csv'), sep=',', header=False)
    return

  def Save_file_heatmap(self, path):
    save_df=self.final_df
    save_path=os.path.join(path,'Heatmapping',f'Heatmap_{self.current_cancer}_{self.current_fpkm}_{self.current_splicing}.csv')
    save_df.to_csv(path_or_buf=save_path, index=False)
# Main code starts here

# For cancers
cancer_name_list=os.listdir('C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/PSI_Data')
cancer_name_list=list(map(lambda x:x[x.rfind('_')+1:x.find('.')],cancer_name_list))
# Previously merged data
fpkm_psi_merge_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Output_FPKM_PSI_match'
# Information about target genes and splicing genes
fpkm_target_gene_list=['fra-1','PPARD','PROM1']
symbol_exon_info_df=pd.read_csv('C:/Users/c/Desktop/LFG_Bioinformatics/Target_list/wnt_symbol_exon_info.txt', sep='\t', header=None)
# File paths for folders containing data.
#fpkm_psi_merged_loc_list=[os.path.join(fpkm_psi_merge_loc, i) for i in os.listdir(fpkm_psi_merge_loc) if i in fpkm_target_gene_list]
# save location
save_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Output_Plot_Data'

main=Separate_top_bottom()

print('Processing Boxplotting')
for fpkm_gene in fpkm_target_gene_list:
  folder_loc=os.path.join(fpkm_psi_merge_loc, fpkm_gene)
  for idx, rows in symbol_exon_info_df.iterrows():
    main.list_of_dfs=[]
    list_row=list(rows)
    splice_gene=list_row[0]
    exon=list_row[1]
    from_exon=float(list_row[2])
    to_exon=float(list_row[3])
    splice_of_interest=f'{splice_gene}_[ES]_"{exon}"({from_exon}->{to_exon})'
    for cancer_name in cancer_name_list:
      file_name=f'{cancer_name}_{fpkm_gene}_{splice_gene}.csv'
      location=os.path.join(folder_loc, file_name)
      main.Separation_boxploting(location=location, fpkm_gene=fpkm_gene, splice_of_interest=splice_of_interest)
    main.Save_file_boxplot(path=save_loc)
'''
print('\nProcessing Heatmapping')
for idx, rows in symbol_exon_info_df.iterrows():
  list_row=list(rows)
  splice_gene=list_row[0]
  exon=list_row[1]
  from_exon=float(list_row[2])
  to_exon=float(list_row[3])
  splice_of_interest=f'{splice_gene}_[ES]_"{exon}"({from_exon}->{to_exon})'
  for cancer_name in cancer_name_list:
    main.current_cancer=cancer_name
    merge_df=main.Separate_heatmapping(location_merge_loc=fpkm_psi_merge_loc,fpkm_gene_list=fpkm_target_gene_list,splice_gene=splice_gene,exon=exon,splice_of_interest=splice_of_interest,cancer_type=cancer_name)
    main.Save_file_heatmap(path=save_loc)
'''

print("Yay! I worked! Debugshrine!")
  