#%%
import os
from datetime import date
import tarfile
import shutil
import psutil
import pandas as pd
import matplotlib as mpl
from pathlib import Path
'''
from lifelines import KaplanMeierFitter
from lifelines.datasets import load_waltons
import gseapy as gp
'''
loc='C:/Users/c/Desktop/LFG_Bioinformatics/Output_Plot_Data/Heatmapping/Heatmapping_Results'
path=Path(__file__).parents[2]
npath=os.path.join(path, 'Output_Plot_Data', 'Heatmapping','Heatmapping_Results')
print(npath)
print(os.listdir(npath))
#%%
