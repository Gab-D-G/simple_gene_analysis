import os
import sys
import pandas as pd
import numpy as np
import nibabel as nb

#csv file with gene_id and gene names
file='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/gene_expression/allen_reference_genes.csv'
gene_info_df=pd.read_csv(file)
info_names=gene_info_df['acronym']
info_id=gene_info_df['id']

target_acronym = sys.argv[1]

if (info_names==target_acronym).sum()==0:
    raise ValueError("The gene acronym wasn't found in the list.")

idx = np.where(info_names==target_acronym)[0][0]
id = str(info_id[idx])


gene_file='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/gene_expression/sagittal_genes.csv'
genes_df=pd.read_csv(gene_file)
roi_id = genes_df['structure_id']

id = str(info_id[idx])
if not id in list(genes_df.columns):
    raise ValueError("The gene is not part of the downloaded list.")

gene_expression = genes_df[str(id)]

atlas_img = nb.load('/data/chamal/projects/Gabriel_DG/gene_analysis_Elisa/resampled_atlas.nii.gz')
data_array = np.asarray(atlas_img.dataobj)
expression_array = np.zeros(data_array.shape)

i=0
multiple_10 = int(len(roi_id)/10)

for roi,expression in zip(roi_id,gene_expression):
    expression_array[data_array==roi] = expression
    if (i%multiple_10==0):
        print(str(int(i/multiple_10)*10)+'% completed.')
    i+=1
nb.Nifti1Image(expression_array, atlas_img.affine, atlas_img.header).to_filename(target_acronym+'_expression.nii.gz')

