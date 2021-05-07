import pandas as pd
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt

map_file='LV1.nii.gz'
map_img=nb.load(map_file)
map_data=map_img.dataobj #gives 3D array with voxel values

atlas_file='resampled_atlas.nii.gz'
atlas_data=nb.load(atlas_file).dataobj

#csv file with ROI ids
label_csv='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/structures.csv'
labels_df=pd.read_csv(label_csv)

data_dict={}
for label in np.asarray(labels_df.get('id'))[1:]: #iterating through the ROI ids
    if np.max(label==atlas_data): #taking a ROI only if it has labeled voxels
        roi_mask=np.asarray(label==atlas_data, dtype=int)  #making a binary mask of the ROI        
        ROI_loading=np.sum(map_data*roi_mask)/np.sum(roi_mask) #taking sum of points on latent variable map within the ROI, then dividing by the number of points in the ROI to get the mean
        data_dict[str(label)]=ROI_loading #assigning the ROI loading to the corresponding ROI id in the dictionary
        print(label)

            
import pickle

with open('ROI_loading.pkl', 'wb') as handle:
    pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

