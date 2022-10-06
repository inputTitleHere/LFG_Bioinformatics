import pandas as pd
from pathlib import Path




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
    canc_row=[]
    #Patient ID
    for names in name_row:
        if '-11A' in names or '-11B' in names:
            continue
        if '-01A' in names or '-01B' in names or '-01C' in names or '-02A' in names:
            canc_row.append(names)
            continue
    #Asterisk * is used to unzip lists.
    sorted_row=[*title_column,*canc_row]
    df=df[sorted_row]
    #df=df.loc[df['Ensembl_ID'].isin(['ENSG00000142871','ENSG00000118523'])]
    df.set_index('Ensembl_ID', inplace=True)
    df=df.loc[['ENSG00000142871','ENSG00000118523']]
    return df, len(title_column), len(canc_row)

def save_to_csv(df, location, index=True):
    df.to_csv(location, sep=',',index=index)
    return


brca_fpkm_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/FPKM_Data/TCGA-BRCA.htseq_fpkm.tsv.gz'
luad_fpkm_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Test_data/FPKM_Data/TCGA-LUAD.htseq_fpkm.tsv.gz'

brca_save_loc='C:/Users/c/Desktop/LFG_Bioinformatics/_output to draw graph/brca_cancut_fpkm.csv'
luad_save_loc='C:/Users/c/Desktop/LFG_Bioinformatics/_output to draw graph/luad_cancut_fpkm.csv'

loc_list=[brca_fpkm_loc,luad_fpkm_loc]
save_loc=[brca_save_loc,luad_save_loc]

for location in loc_list:
    cancer_name=location[location.rfind('-')+1:location.find('.')]
    save_path=[i for i in save_loc if cancer_name.lower() in i]
    df,*_=open_fpkm(location)
    save_to_csv(df, save_path[0])