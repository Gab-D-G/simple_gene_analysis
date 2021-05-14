import pandas as pd
import numpy as np
import nibabel as nb
import pickle

with open('ROI_loading.pkl', 'rb') as handle:
    roi_dict=pickle.load(handle)

gene_file='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/gene_expression/sagittal_genes.csv'
genes_df=pd.read_csv(gene_file)
genes_array=np.asarray(genes_df) #array format of the gene expression data to extract rows

data_header=list(genes_df.columns)
data_header.append('ROI_loading')
data_df=pd.DataFrame(columns=data_header) # create a new dataframe to store the ROI loading on latent variable map, for each associated roi id

for roi_id in list(roi_dict.keys()): #will iterate through the different roi with their associated gene expression
    roi_loading=roi_dict[roi_id] #get the latent variable loading for the roi
    try:
        roi_idx=np.where(genes_df['structure_id']==int(roi_id))[0][0]
    except:
        continue
    roi_data=list(genes_array[roi_idx, :]) #make a list with gene expression values for that given roi
    roi_data.append(roi_loading) #add the roi loading to the list
    array=np.empty([1, len(data_header)]) #need to transform the data to a rowxcolumn array to create a DataFrame
    array[0,:]=np.asarray(roi_data) #insert the data for that roi into the array
    nan_idx=np.where(np.isnan(array)) #get index where there is no gene expression values (with 'nan' values)
    array[nan_idx]=0 #set these values to 0, since linear regression probably won't work with 'nan' values
    roi_df=pd.DataFrame(array, index=[1], columns=data_header) #create a data frame with a single row, representing gene expression/latent variable loading data for the given ROI
    data_df=data_df.append(roi_df) #append the roi info to the overall DataFrame

'''
Evaluate p-values for all genes with univariate linear regression
'''
import statsmodels.api as sm

stat_dict={}
gene_columns=data_df.columns[1:-1]  #skip the first column which corresponds to structure ID, and don't select last column which is the ROI loading
for gene in gene_columns:
    X = data_df[gene]
    Y = data_df['ROI_loading']
    #X = sm.add_constant(X) # can add an intercept to the linear regression if desired
    model = sm.OLS(Y, X).fit() #fitting a linear regression with Ordinary Least Square, where gene expression is a predictor of ROI loading on the latent variable map
    stat_dict[gene]=[model.tvalues[0], model.pvalues[0]]


#csv file with gene_id and gene names
file='/data/chamal/projects/Gabriel_DG/atlases/AllenBrain/gene_expression/allen_reference_genes.csv'
gene_info_df=pd.read_csv(file)
info_names=gene_info_df['acronym']
info_id=gene_info_df['id']


t_values=[]
p_values=[]
gene_names=[]
keys=stat_dict.keys()
for i in range(len(info_id)):
    id=str(info_id[i]) #converting the integer id to string since the keys of stat_dict are strings
    if id in keys:
        t_values.append(stat_dict[id][0])
        p_values.append(stat_dict[id][1])
        gene_names.append(info_names[i])

#convert all list to an array to be able to use np.argsort
gene_names_array=np.asarray(gene_names)
t_values_array=np.asarray(t_values)
p_values_array=np.asarray(p_values)

#get the index to sort the gene_id and p-values arrays in increasing order of the associated t-values
sort_index=np.argsort(t_values_array)
sorted_t_values=t_values_array[sort_index]
sorted_p_values=p_values_array[sort_index]
sorted_gene_names=gene_names_array[sort_index]

data=np.asarray([sorted_gene_names, sorted_p_values]).transpose() #convert the two vector to a single array, and transpose to proper dimensions to create a Dataframe

#Create Dataframe with one column for the gene name, and one column for the associated p-value, where the first values are the most significant negative associations and last are the most significant positive associations
df = pd.DataFrame(data, columns = ['Gene Name', 'p-values'])
df.to_csv('LV_gene_associations.csv')
