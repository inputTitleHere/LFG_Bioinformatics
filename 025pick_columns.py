import os
from numpy.lib.npyio import save 
import pandas as pd
from pathlib import Path

gene_name='APC'
path=os.path.join(Path(__file__).parents[1],'Test_Output_2',f'{gene_name}')

za_list=os.listdir(path)
for files in za_list:
    loc=os.path.join(path, files)
    df=pd.read_csv(loc, sep=',')
    n_df=df[['symbol','exons','mean_diff','pos_neg']]
    save_path=os.path.join(path, f'le_{files}')
    n_df.to_csv(save_path,sep=',',index=False)