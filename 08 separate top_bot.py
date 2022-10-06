from typing import final
from numpy.core.numeric import full
import pandas as pd
import os
from pathlib import Path

from pandas.core.reshape.merge import merge

# Read files from test_output_4 to separate top and bottom 20% and save it in separate csv files

def open_csv_file(path, gene_symbol, cancer_name, splice_gene, exon, from_exon, to_exon, percentage=20):
    '''path : Set path to test_ouput_4 \ngene_symbol : Name of pathway target gene\ncancer_name : Put target cancer name here \nexon, from_exon, to_exon : Put exon information as numbers \npercentage : Sets cutoff location. Default is 20'''
    print(cancer_name)
    file_path=os.path.join(path, gene_symbol)
    file_list=os.listdir(file_path)
    file_name=[n for n in file_list if f'{cancer_name}_{gene_symbol}_{splice_gene}' in n]
    if len(file_name)==1:
        file_loc=os.path.join(file_path, file_name[0])
        exon_info=f'"{exon}"({from_exon}->{to_exon})'
    #===#
    if percentage>1:
        percentage=percentage*0.01
    #===#
    main_df=pd.read_csv(file_loc, sep=',', encoding='utf-8')
    index_list=list(main_df.columns)
    try:
        #Check if list is empty
        splice_of_interest=[i for i in index_list if exon_info in i][0]
    except IndexError:
        print(f'{cancer_name} has no {splice_gene}_{exon_info}')
        return
    main_df=main_df[splice_of_interest].copy()
    locus=round(main_df.size*percentage)
    top_side_df=main_df.iloc[0:locus]
    bot_side_df=main_df.iloc[main_df.size-locus:main_df.size]
    #====#
    # Renaming #
    #====#
    top_side_df=top_side_df.to_frame()
    bot_side_df=bot_side_df.to_frame()
    full_exon_info=top_side_df.columns[0]
    top_side_df.rename(columns={full_exon_info:f'{cancer_name}_{gene_symbol}_top_{percentage*100}%'},inplace=True)
    bot_side_df.rename(columns={full_exon_info:f'{cancer_name}_{gene_symbol}_bot_{percentage*100}%'}, inplace=True)
    top_side_df.reset_index(drop=True,inplace=True)
    bot_side_df.reset_index(drop=True,inplace=True)
    merge_df=pd.concat([top_side_df, bot_side_df],axis=1)
    return merge_df

def save_df(path, gene_symbol, cancer_name, splice_gene, exon_info):    
    return

#=======================#
# Main code starts here #
#=======================#

folder='C:/Users/c/Desktop/LFG_Bioinformatics/Output_FPKM_PSI_match'
symbol_exon_info_df=pd.read_csv('C:/Users/c/Desktop/LFG_Bioinformatics/Target_list/wnt_symbol_exon_info.txt', sep='\t')
gene_symbols=['fra-1','PROM1','PPARD']
cancer_name_list=os.listdir('C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/FPKM_Data')
cancer_name_list=list(map(lambda x:x[x.find('-')+1:x.find('.')],cancer_name_list))
print(cancer_name_list)

for idx,rows in symbol_exon_info_df.iterrows():
    list_row=list(rows)

    splice_gene=list_row[0]
    exon=list_row[1]
    from_exon=float(list_row[2])
    to_exon=float(list_row[3])
    for gene_symbol in gene_symbols:
        if ':' in exon:
            full_exon_name=f'{gene_symbol}={splice_gene}_{exon}({from_exon}-{to_exon})'
            full_exon_name=full_exon_name.replace(':','+')
        else:
            full_exon_name=f'{gene_symbol}={splice_gene}_{exon}({from_exon}-{to_exon})'

        print(full_exon_name)
        merge_dfs=[]
        for cancer_name in cancer_name_list:
            merge_dfs.append(open_csv_file(folder, gene_symbol=gene_symbol, cancer_name=cancer_name, splice_gene=splice_gene,exon=exon, from_exon=from_exon, to_exon=to_exon))
        final_df=pd.concat(merge_dfs, axis=1,ignore_index=False)
        print(final_df)
        final_df.to_csv(path_or_buf=f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_5/{full_exon_name}.csv',index=False,encoding='utf-8')