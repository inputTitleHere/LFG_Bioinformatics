from numpy.lib.npyio import save
import pandas as pd
import seaborn as sns
import os
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from seaborn.external.husl import f

# By using defs, this code aims to draw a heatmap once at a time.

class Heatmapping:
  def __init__(self):
    self.df=pd.DataFrame()
    self.heatmap_df_list=[]
    self.colorset=None
    self.readlocation=None
    self.filename=None
    self.savelocation=None

  def open_dataframe(self, path):
    '''Reads comma (',') separated csv and returns a dataframe'''
    self.df=pd.read_csv(path, sep=',', encoding='utf-8',index_col=False)
    self.readlocation=path
    self.filename=path[path.rfind('\\')+1:path.find('.')]
    self.savelocation=path[:path.find('\\')]
    return

  def sep_dataframe(self,fpkm_gene_list):
    '''Get proper data that will be used to draw a heatmap'''
    column_names=list(self.df.columns)
    #Add first column for splicing gene. See the file structure in 'Output_Plot_data/Heatmapping' folder
    self.heatmap_df_list.append(self.df[column_names[0]])
    for headername in column_names:
      if headername in fpkm_gene_list:
        self.heatmap_df_list.append(self.df[headername])
    return

  def draw_heatmap(self):
    '''Draws heatmap from loaded dataframe. Current version only supports three column heatmap'''
    # Set colormap of Red - Black - Green
    color=[(1,0,0),(0,0,0),(0,1,0)]
    color_map=LinearSegmentedColormap.from_list(name='red_black_green',colors=color,N=300)
    # Three figures
    fig, (ax1,ax2,ax3)=plt.subplots(ncols=3, figsize=(6,10))
    fig.subplots_adjust(wspace=0)
    splice_df, second_df, third_df=self.heatmap_df_list
    splice_df=splice_df.to_frame()
    second_df=second_df.to_frame()
    third_df=third_df.to_frame()
    #=== Plotting ===#
    sns.heatmap(splice_df,cmap=color_map, ax=ax1,yticklabels=False,cbar=False)
    cbar1=fig.colorbar(ax1.collections[0],ax=ax1,orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    ax1.set_ylabel('')
    #===Second df HEATMAP===#
    #Manpulate vmin and vmax as boundaries for the heatmap plot
    sns.heatmap(second_df,vmin=0, vmax=4, cmap=color_map,ax=ax2, yticklabels=False, cbar=False)
    cbar2=fig.colorbar(ax2.collections[0],ax=ax2, orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    cbar2.set_ticks([n for n in range(0,8)])
    cbar2.set_ticklabels([n for n in range(0,8)])
    ax2.set_ylabel('')
    #===Third df HEATMAP===#
    sns.heatmap(third_df,vmin=2, vmax=6, cmap=color_map,ax=ax3, yticklabels=False, cbar=False)
    cbar3=fig.colorbar(ax3.collections[0],ax=ax3, orientation='horizontal',use_gridspec=False, pad=0.06,shrink=0.9,aspect=5)
    labels3=[n for n in range(2,7)]
    cbar3.set_ticks(labels3)
    cbar3.ax.set_xticklabels(labels3, rotation=0)
    ax3.set_ylabel('')
    return
  
  def save_to_svg(self):
    loc=self.savelocation
    save_loc=os.path.join(loc, 'Heatmapping_Results')
    if not os.path.exists(save_loc):
      os.makedirs(save_loc)
    save_path=os.path.join(save_loc, str(self.filename)+'.svg')
    plt.savefig(save_path, format='svg')
    return

#================= Main ===================#

if __name__=='__main__':
  heatmap_data_loc='C:/Users/c/Desktop/LFG_Bioinformatics/Output_Plot_Data/Heatmapping'
  heatmap_data_list=os.listdir(heatmap_data_loc)
  #=== Parameters ===#
  cancer_type=['BRCA']
  splice_gene=['RPS24','ERBB2IP']
  fpkm_gene_list=['fra-1','PPARD']
  #=== Parameters ===#
  target_df_list=[i for i in heatmap_data_list if i[i.find('_')+1:i.find('__')] in cancer_type]
  target_df_list=[i for i in target_df_list if i[i.find('__')+2:i.find('.')] in splice_gene]
  for csv in target_df_list:
    H=Heatmapping()
    H.open_dataframe(os.path.join(heatmap_data_loc, csv))
    H.sep_dataframe(fpkm_gene_list=fpkm_gene_list)
    H.draw_heatmap()
    H.save_to_svg()

  