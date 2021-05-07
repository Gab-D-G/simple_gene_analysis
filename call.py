# get gene and structure IDs; if you want to provide structure IDs you can do so, but
# by default I'd suggest the Rubinov et al., 2015, PNAS set
import os
from abagen import mouse
import pandas as pd
import pickle
genes = mouse.fetch_allenref_genes('id')
#structures = mouse.fetch_rubinov2015_structures('id')

import pandas as pd
import numpy as np

atlas='/data/chamal/projects/Gabriel_DG/dFC_project/Grandjean_2014/pip_output_20190401/sub-mi272871_ses-1_run-1_anat_labels.nii.gz'
atlas='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/ABI_labels_on_DSURQE_40micron.nii.gz'

import nibabel as nb
atlas_data=nb.load(atlas).dataobj
atlas_img=nb.load(atlas)

label_csv='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/structures.csv'
labels_df=pd.read_csv(label_csv)

label_list=[]
for label, hierarchy in zip(np.asarray(labels_df.get('id'))[1:], np.asarray(labels_df.get('structure_id_path'))[1:]):
    #if np.max(label==atlas_data): #taking a ROI only if it has labeled voxels
    if '/997/' in hierarchy: #download for all regions in the adult mouse atlas
        label_list.append(label)

structures=label_list

# get expression data for provided genes; this will take a while since there are
# ~20000 genes. we'll query ~100 genes at a time (querying all at once will error)

expression=[]
iter=0

chunk = 10
for i in range(iter, len(genes), chunk):
    try:
        expression += [mouse.get_unionization_from_gene(id=genes[i:i + chunk],
                                                        structures=structures,
                                                        slicing_direction='sagittal')] #taking only sagittal data
    except:
        print("No results for genes "+str(genes[i:i + chunk]))

    if i%1000==0:
        print(i)
        with open('temp_genes.pkl', 'wb') as handle:
            pickle.dump({'expression':expression,'iter':i+chunk}, handle, protocol=pickle.HIGHEST_PROTOCOL)


# now concatenate across all genes and pivot so that it's a structure x gene dataframe
expression_df = pd.pivot_table(pd.concat(expression).reset_index(),
                             index='structure_id', columns='gene_id',
                             values='expression_density')

with open('sagittal_genes.pkl', 'wb') as handle:
    pickle.dump(expression_df, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('sagittal_genes.pkl', 'rb') as handle:
    expression_df=pickle.load(handle)

expression_df.to_csv('sagittal_genes.csv')

