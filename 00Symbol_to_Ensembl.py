import mygene
from pathlib import Path
import os


mg=mygene.MyGeneInfo()

#xli = ['CTNNB1','APC','Axin1','NkD1','YAP1','TAZ']
#gene_list=['GSK3b','DVL1','LEF1']
file_path=Path(__file__).parents[1]
file_loc=os.path.join(file_path,'Target_list')
save_loc=os.path.join(file_path,'Target_to_Ensembl')

#input target keyword here
keyword='hippo'
target_file_name=''
file_name=''
for files in os.listdir(file_loc):
    if keyword.lower() in files.lower():
        target_file_name=files
        
        
        
        
        
genelist=[]
with open(os.path.join(file_loc, target_file_name), mode='r', encoding='utf-8') as t:
    gene_list=t.readlines()

using_gene_list=gene_list

result=mg.querymany(using_gene_list, scopes='symbol', fields='ensembl.gene', species='human')
print_list=[]
for genes in result:
    printing_string=genes['query']+'\t'+genes['ensembl']['gene']
    print_list.append(printing_string)
    print(printing_string)

#save here

with open(os.path.join(save_loc, f'{target_file_name[:target_file_name.rfind(".")]}_Ensembl.txt'), mode='w+', encoding='utf-8') as save_txt:
    save_txt.writelines(f"{line}\n" for line in print_list)

print('Finished!')

