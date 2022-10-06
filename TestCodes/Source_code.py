import pandas as pd
import numpy as np
import zipfile
import tarfile
import shutil
import psutil
import gzip
import os
from pathlib import Path
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter as kmf
from datetime import date
from statistics import mean
'''
Todo List = [1. Using shutil.move(), make code that moves processed files to another location. preventing accidental clashes among other samples.]
'''
#====================================================================#
#                           Coding for FPKM                          #
#====================================================================#
def get_FPKM_file_location():
    root_path=Path(__file__).parents[1]
    folders=os.listdir(root_path)
    for folder in folders:
        if 'fpkm' in folder.lower() and 'input' in folder.lower():
            folder_path=folder
    path_FPKM=os.path.join(root_path,folder_path)
    file_lst=os.listdir(path_FPKM)
    file_loc_lst=[]
    for file_name in file_lst:
        if file_name.startswith('~$'):
            continue
        if file_name.endswith(".tsv") or file_name.endswith('.gz'):
            file_loc=os.path.join(path_FPKM, file_name)
            file_loc_lst.append(file_loc)
    return file_loc_lst
# return Example:['c:\\Users\\c\\Desktop\\LFG_Bioinformatics\\input_data_FPKM\\TCGA-LUAD.htseq_fpkm.tsv', 'c:\\Users\\c\\Desktop\\LFG_Bioinformatics\\input_data_FPKM\\TCGA-LUAD.htseq_fpkm.tsv.gz']

def open_FPKM(file_loc):
    '''This method takes single locations. If a list is returned by get_FPKM_file_location(), then use a loop.'''
    if file_loc.endswith('.gz'):
        df=pd.read_csv(file_loc, sep='\t',compression='gzip', index_col=False)
        file_name=os.path.basename(file_loc)
        return df, file_name
    elif file_loc.endswith('.tsv'):
        df=pd.read_csv(file_loc, sep='\t', index_col=False)
        file_name=os.path.basename(file_loc)
        return df, file_name

def process_dataframe(dataframe, pvalue_limit=0.05):
    #get first row
    df=dataframe
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
    title_count=len(title_column)
    norm_count=len(norm_row)
    canc_count=len(canc_row)
    df=df[sorted_row]
    # calculate pvalue and mean value from here.
    pval_lst=[]
    norm_mean_lst=[]
    canc_mean_lst=[]
    fold_lst=[]
    norm_count_lst=[]
    canc_count_lst=[]
    name_lst=[]
    #This is where it take a lot of time.
    #for loops are sloooooooooow
    for _, series in df.iterrows():
        norm_lst=series[title_count:title_count+norm_count]
        canc_lst=series[title_count+norm_count:]
        norm_lst=((2**norm_lst)-1)
        norm_lst=norm_lst.loc[(norm_lst!=0)]
        canc_lst=((2**canc_lst)-1)
        canc_lst=canc_lst.loc[(canc_lst!=0)]
        norm_total_count=len(norm_lst)
        canc_total_count=len(canc_lst)
        if norm_total_count<norm_count*0.5 or canc_total_count<canc_count*0.5:
            continue
        norm_mean=mean(norm_lst)
        canc_mean=mean(canc_lst)
        if norm_mean<1 or canc_mean<1:
            continue
        result=ttest_ind(norm_lst,canc_lst, equal_var=False)
        pval=result.pvalue
        if pval>pvalue_limit:
            continue
        # Append values to list and convert it to a new dataframe
        # Series are [[]] shaped. Get EnsemblID.
        ensembl_symb=list(series[:title_count])[0]
        try:
            name_string=ensembl_symb[:ensembl_symb.index('.')]
            name_lst.append(name_string)
        except ValueError:
            pass
        pval_lst.append(pval)
        norm_mean_lst.append(norm_mean)
        canc_mean_lst.append(canc_mean)
        fold=canc_mean/norm_mean
        fold_lst.append(fold)
        norm_count_lst.append(norm_total_count)
        canc_count_lst.append(canc_total_count)
    dic={'Ensembl_ID':name_lst,'p_value':pval_lst,'Fold':fold_lst,'Normal_mean':norm_mean_lst,'Cancer_mean':canc_mean_lst,'Normal_count':norm_count_lst,'Cancer_count':canc_count_lst}
    calculated_df=pd.DataFrame(dic)
    calculated_df=calculated_df.sort_values(by='Fold',axis=0, ascending=False)
    return calculated_df
#==================================#
def get_clinical_list(dataframe, percentage=25, ensembl_id_list=None):
    '''Dataframe is for FPKM.'''
    dictionary={}
    if ensembl_id_list is not None:
        df=dataframe[ensembl_id_list]
        # If using this part, might want to deal with multiple ensemble ids per gene retrieved from mygene.
    else:
        df=dataframe
    df.set_index('Ensembl_ID',drop=True,inplace=True)
    columns_count=len(df.columns)
    portion=round(columns_count*(percentage/100))
    for index,series in df.iterrows():
        if series.mean() < 1:
            continue
        temp_series=series[series.values!=0]
        temp_sorted=series.sort_values()
        temp_column_name_lst=temp_series.index.tolist()
        if len(temp_column_name_lst)<portion:
            continue
        dictionary[series.name]=temp_column_name_lst
    return dictionary, portion
    
        

        
        

    print(df)
    print(columns_count)

def save_to_Excel(dataframe, file_name, file_type='csv'):
    root_path=Path(__file__).parents[1]
    new_file_name=f'{date.strftime(date.today(),"%m%d%y")}_{file_name[:file_name.index(".tsv")]}'
    save_folder=os.path.join(root_path,'Output_Data')
    save_loc=os.path.join(save_folder, new_file_name)
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    if file_type.lower()=='excel' or file_type.lower()=='e':
        save_loc=f'{save_loc}.xlsx'
        with pd.ExcelWriter(save_loc, mode='w', engine='openpyxl') as writer:
            dataframe.to_excel(writer, sheet_name='Sheet1', index=False)
    if file_type.lower()=='csv' or file_type.lower()=='c':
        save_loc=f'{save_loc}.csv'
        dataframe.to_csv(save_loc, sep=',',index=False)

    print(f"File '{new_file_name}'' is saved at '{save_folder}'.")
    return

def remove_normal(dataframe):
    name_row=list(dataframe.columns)
    title_column=[name_row[0]]
    canc_row=[]
    for names in name_row:
        if '-01A' in names or '-01B' in names or '-01C' in names or '-02A' in names:
            canc_row.append(names)
    sorted_row=[*title_column, *canc_row]
    # 명시적 데이터프레임을 copy() 를 이용하여 제작한다.  SettingWithCopyWarning을 방지.
    df=dataframe[sorted_row].copy()
    new_ensembl_frame=df.loc[:,'Ensembl_ID'].map(lambda x : x[:x.index('.')])
    df['Ensembl_ID']=new_ensembl_frame
    df=df.rename(columns={'Ensembl_ID':'Ensembl_ID-'})
    df=df.rename(columns=lambda x : x[:x.rfind('-')])
    return df
#====================================================================#
#                           Coding for PSI                          #
#====================================================================#

def get_PSI_file_location():
    root_path=Path(__file__).parents[1]
    folders=os.listdir(root_path)
    for folder in folders:
        if 'psi' in folder.lower() and 'input' in folder.lower():
            folder_path=folder
    path_PSI=os.path.join(root_path, folder_path)
    file_lst=os.listdir(path_PSI)
    file_loc_lst=[]
    for file_name in file_lst:
        if file_name.endswith('.zip') or file_name.endswith('.txt'):
            file_loc=os.path.join(path_PSI,file_name)
            file_loc_lst.append(file_loc)
    print(file_loc_lst)
    return file_loc_lst



#====================================================================#
#                           Coding for Clinical                          #
#====================================================================#


def get_Clinical_file_location():
    root_path=Path(__file__).parents[1]
    folders=os.listdir(root_path)
    for folder in folders:
        if 'clinical' in folder.lower() and 'input' in folder.lower():
            folder_path=folder
    path_Clinical=os.path.join(root_path, folder_path)
    file_lst=os.listdir(path_Clinical)
    file_loc_lst=[]
    for file_name in file_lst:
        if file_name.endswith('.gz') and 'clinical' in file_name.lower():
            file_loc=os.path.join(path_Clinical, file_name)
            file_loc_lst.append(file_loc)
        if file_name.endswith('.tsv') and 'clinical' in file_name.lower():
            file_loc=os.path.join(path_Clinical, file_name)
            file_loc_lst.append(file_loc)
    return file_loc_lst

def open_Clinical(file_loc):
    if file_loc.endswith('.gz') and 'clinical' in file_loc.lower():
        with tarfile.open(file_loc) as tf:
            members=tf.getnames()
            for files in members:
                if 'clinical' in files:
                    clinical_file=tf.extractfile(files)
            df=pd.read_csv(clinical_file, sep='\t',index_col=False)
        file_name=os.path.basename(file_loc)
        return df, file_name

def remove_Clinical_Duplicates(dataframe):
    df=dataframe
    df=df.drop_duplicates(subset='case_submitter_id')
    df.reset_index(inplace=True, drop=True)
    return df

def extract_Clinical(dataframe):
    column_lst=['case_submitter_id','days_to_death','vital_status','ajcc_pathologic_stage','days_to_last_follow_up']
    df=dataframe[column_lst].copy()
    return df

def match_Clinical_to_Others(clinical_df, fpkm_df=None,psi_df=None):
    '''The centerpiece of clinical process'''
    if fpkm_df is not None:
        # We have both clinical and FPKM dataframes.
        fpkm_df=remove_normal(fpkm_df)
        patient_by_ensembl_dic, int_portion=get_clinical_list(fpkm_df,percentage=25)
        c_df=extract_Clinical(remove_Clinical_Duplicates(clinical_df))
        new_clinical_df=move_Clinical_Time(c_df)
        

        return c_df
    if psi_df is not None:
        print('running psi_df')
        pass

def clinical_top_bot_divider(ensembl_dic, int_portion):
    '''Divides top portion and bottom portion into separate dictionaires.'''
    top_dic={}
    bot_dic={}
    for ensembl_ID in ensembl_dic:
        patient_lst=ensembl_dic[ensembl_ID]
        top=patient_lst[:int_portion]
        bot=patient_lst[len(patient_lst)-int_portion:]
        top_dic[ensembl_ID]=top
        bot_dic[ensembl_ID]=bot
    return top_dic,bot_dic


def move_Clinical_Time(dataframe):
    df=dataframe
    row=[]
    for _,series in df.iterrows():
        if series['vital_status']=='Alive':
            row.append(series['days_to_last_follow_up'])
        elif series['vital_status']=='Dead':
            row.append(series['days_to_death'])
    df['merged_days']=pd.Series(row)
    df=df.loc[df['merged_days']!="'--"]
    return df

def save_Clinical(dataframe, file_name, file_type='csv'):
    root_path=Path(__file__).parents[1]
    new_file_name=f'{date.strftime(date.today(),"%m%d%y")}_{file_name[:file_name.rindex(".")]}'
    save_folder=os.path.join(root_path,'Output_Data')
    save_loc=os.path.join(save_folder, new_file_name)
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    if file_type.lower()=='excel' or file_type.lower()=='e':
        save_loc=f'{save_loc}.xlsx'
        with pd.ExcelWriter(save_loc, mode='w', engine='openpyxl') as writer:
            dataframe.to_excel(writer, sheet_name='Sheet1', index=False)
    if file_type.lower()=='csv' or file_type.lower()=='c':
        save_loc=f'{save_loc}.csv'
        dataframe.to_csv(save_loc, sep=',',header=True,index=False )

    print(f"File '{new_file_name}'' is saved at '{save_folder}'.")
    return
        





    






#====================================================================#
#                           Coding for Main run .                          #
#====================================================================#




fpkm_run=True
psi_run=False
clinical_run=False

if fpkm_run:
    print('======Running FPKM File=====')
    file_lst=get_FPKM_file_location()
    for file_locs in file_lst:
        dataframe, file_name=open_FPKM(file_locs)
        dataframe=process_dataframe(dataframe=dataframe)
        save_to_Excel(dataframe=dataframe, file_name=file_name)

if psi_run:
    print('======Running PSI File=====')
    get_PSI_file_location()

if clinical_run:
    print('======Running Clinical File=====')
    clinical_file_lst=get_Clinical_file_location()
    fpkm_file_loc_lst=get_FPKM_file_location()
    fpkm_df, _=open_FPKM(fpkm_file_loc_lst[0])
    for file_loc in clinical_file_lst :
        #Just put one in.
        clinical_df, clinical_file_name=open_Clinical(file_loc)
        c_df=match_Clinical_to_Others(clinical_df, fpkm_df=fpkm_df)
        #print(c_df)
        #save_Clinical(c_df,clinical_file_name)



        print('memory usage =', psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)



print(f'[Execution Report : FPKM_RUN={fpkm_run}\tPSI_RUN={psi_run}\tClinical_RUN={clinical_run}]')