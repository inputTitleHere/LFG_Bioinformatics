#%%
import os
from matplotlib.colors import LinearSegmentedColormap
from pandas.core.accessor import register_index_accessor 
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


#===#
cancer_name_list=['LUAD','BRCA','STAD','BLCA','OV','PRAD']

#===#
'''
brca_ccn1_path='C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN1/LUAD_CCN1.csv'
brca_ccn2_path='C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN2/LUAD_CCN2.csv'

brca_ccn1_df=pd.read_csv(brca_ccn1_path, sep=',', encoding='utf-8')
brca_ccn2_df=pd.read_csv(brca_ccn2_path, sep=',', encoding='utf-8')
'''
#===#
exons, from_e, to_e=17,float(16),float(18)
splice_info=f'"{exons}"({from_e}->{to_e})'
#===#
for names in cancer_name_list:
    print(f'___Processing {names}___\n\n\n')
    file_path_ccn1=f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN1/{names}_CCN1.csv'
    file_path_ccn2=f'C:/Users/c/Desktop/LFG_Bioinformatics/Test_Output_4/CCN2/{names}_CCN2.csv'
    df_ccn1=pd.read_csv(file_path_ccn1, sep=',', encoding='utf-8')
    df_ccn2=pd.read_csv(file_path_ccn2, sep=',', encoding='utf-8')
    splice_type=[i for i in df_ccn1.columns if splice_info in i]
    df_ccn1=df_ccn1[~df_ccn1[f'{splice_type[0]}'].isnull()]
    df_ccn2=df_ccn2[~df_ccn2[f'{splice_type[0]}'].isnull()]
    df_ccn1.sort_values(splice_type[0], inplace=True, ascending=False)
    df_ccn2.sort_values(splice_type[0], inplace=True, ascending=False)
    df_ccn1.set_index('patient_id', drop=True, inplace=True)
    df_ccn2.set_index('patient_id', drop=True, inplace=True)
    indexes=list(df_ccn1.index)
    df_ccn2=df_ccn2.loc[indexes]

    #=========================================================================#
    ccn1_df=df_ccn1['CCN1']
    ccn1_df=ccn1_df.to_frame()
    ccn2_df=df_ccn2['CCN2']
    ccn2_df=ccn2_df.to_frame()
    #===#
    splice_columns=df_ccn1.columns
    splice_of_interest=[i for i in splice_columns if splice_info in i]
    splice_df=df_ccn1[splice_of_interest[0]]
    splice_df=splice_df.to_frame()
    #=========================================================================#
    color_plate=sns.color_palette('vlag',as_cmap=True)

    color=[(1,0,0),(0,0,0),(0,1,0)]
    color_map=LinearSegmentedColormap.from_list(name='red_black_green',colors=color,N=300)

    fig, (ax1,ax2,ax3)=plt.subplots(ncols=3, figsize=(6,10))
    fig.subplots_adjust(wspace=0)
    fig.suptitle(f'{names}', fontsize=15)
    #===#
    sns.heatmap(splice_df,cmap=color_map, ax=ax1,yticklabels=False,cbar=False)
    cbar1=fig.colorbar(ax1.collections[0],ax=ax1,orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    ax1.set_ylabel('')
    #===CCN1 HEATMAP===#
    sns.heatmap(ccn1_df,vmin=2, vmax=10, cmap=color_map,ax=ax2, yticklabels=False, cbar=False)
    cbar2=fig.colorbar(ax2.collections[0],ax=ax2, orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    cbar2.set_ticks([n for n in range(2,11)])
    cbar2.set_ticklabels([n for n in range(2,11)])
    ax2.set_ylabel('')
    #===CCN2 HEATMAP===#
    sns.heatmap(ccn2_df,vmin=2, vmax=10, cmap=color_map,ax=ax3, yticklabels=False, cbar=False)
    cbar3=fig.colorbar(ax3.collections[0],ax=ax3, orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    labels3=[n for n in range(2,11)]
    cbar3.set_ticks(labels3)
    cbar3.ax.set_xticklabels(labels3, rotation=0)
    ax3.set_ylabel('')
    #===#

    save_path=os.path.join(Path(__file__).parents[1], 'Plot_results',f'{names}_CCN1_CCN2_MARK3_{exons}.svg')
    plt.savefig(save_path, format='svg')
    plt.show()
#graph=sns.heatmap(data=brca_ccn1_df,cmap=sns.color_palette("Spectral",as_cmap=True),yticklabels=False)
#graph.set(ylabel=None)
#

# %%
