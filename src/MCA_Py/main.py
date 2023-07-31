import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy, mca, dist

"""
Ramin Mohammadi. main.py: File showing how to use functions in MCA package. Shows example of reading data from csv files and storing them in a pandas data frame, calling RunMCA() to perform MCA (Multiple Correspondence Analysis), and calling GetDistances() to get distance between resulting  coordinates from MCA.
    Copyright (C) 2023  Ramin Mohammadi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses """




# must pip install scanpy, numpy, pandas, matplotlib


##### if wanting to use scaled data (from .h5ad file) ##### -> Refer to loadingInputData.txt to get and use the scaled data from a .rds dataset stored in a seurat object, BUT this data does not seem to be used when running in R (rather the assay of the data is used) so IGNORE THIS COMMENTED STEP
# # In dataset, ROWS should be genes and COLUMNS should be cells.
# annDataObject = scanpy.read_h5ad("BaronScaledNorm.h5ad") # filename will differ but should be .h5ad extension
# annDataObject = annDataObject.transpose() # transpose to get rows to be genes and columns to be cells
# # turn AnnData object into pandas data frame
# scaled_data_df = pd.DataFrame(annDataObject.X, index=annDataObject.obs_names, columns=annDataObject.var_names)
#########################################


#IMPORTANT NOTE for below code: Refer to loadAssay.txt to see how to get the required .csv files and by using the data in a .rds dataset that contains a seurat object
# (you will locally have to create the .csv files yourself using the step described in loadAssay.txt due to the large size of the datasets and the .csv files generated)

# Assay of dataset
assay_df = pd.read_csv('Baron_assay_data.csv') 
# row names of assay
row_names = pd.read_csv('Baron_assay_rownames.csv')

# set row names to data frame's row names
assay_df.index = np.ndarray.flatten(row_names.values)

# If want to test with small input, use assay_df_smaller 
# (change the index values as needed). Currently uses first 400 rows and first 200 columns
# assay_df_smaller = assay_df.iloc[:400, :200]

# Perform MCA on input data
# result is an object containing:
    # featuresCoordinates: pandas data frame where rows are genes and columns are nmcs (default value for nmcs is 50) different dimensions
    # cellsCoordinates: pandas data frame where rows are cells and columns are nmcs (default value for nmcs is 50) different dimensions
    # stdev: numpy array containing the singular values created during MCA when performing Singular Value Decomposition
result = mca.RunMCA(assay_df)

# IMPORTANT NOTE: if specifying a value for nmcs in RunMCA(), must provide same value for dims when calling GetDistances()
# Run MCA using only specified features (here is first 60 features) and specifying nmcs
#result = mca.RunMCA(assay_df, features=assay_df.index.tolist()[:60], nmcs=30)

# MCA on smaller input data
# result = mca.RunMCA(assay_df_smaller)

# Optionally print results
# print("Singular values (stdev): ", result.stdev)
# print("Features coordinates: ", result.featuresCoordinates.shape, result.featuresCoordinates)
# print("Cells coordinates: ", result.cellsCoordinates.shape, result.cellsCoordinates)


# Optionally write result coordinates in form of data frames to CSV files to see entire results
# result.featuresCoordinates.to_csv('Results_csv/featuresCoordinates.csv', index=False)
# result.cellsCoordinates.to_csv('Results_csv/cellsCoordinates.csv', index=False)

# Optionally plot coordinates of first 2 dimensions (use first 2 columns in each coordinate data frame as x and y) 
# IMPORTANT NOTE: 
    # The orientation of the coordinates plotted is not consistent meaning the absolute values of the coordinates, if ran multiple times, will be the same, but the signs of the coordinate values may differ meaning one time the value will be positive and another time will be negative.
    # BUT, with this in mind, the shape of the scatterplot, no matter the orientaton, is always the same meaning the distance between gene X and cell Y is still always the SAME. Thus, does not affect the next step in cell identification being to find the distance between genes and cells in an nmcs dimensional space.
# plt.scatter(result.featuresCoordinates.iloc[:, 0], result.featuresCoordinates.iloc[:, 1], label='Genes(features) Coordinates', c='blue', marker='x')
# plt.scatter(result.cellsCoordinates.iloc[:, 0], result.cellsCoordinates.iloc[:, 1], label='Cells Coordinates', c='red', marker='o', alpha=0.5)
# # Add labels and legend
# plt.xlabel('MCA_1')
# plt.ylabel('MCA_2')
# plt.title('MCA Scatter Plot')
# plt.legend()
# # Show the plot
# plt.show()

# Find DISTANCES between the coordinates of features (genes) and cells
# CellGeneDistances will be a pandas data frame, 2d matrix containing distances between genes and cells (rows are genes, cells are columns)
CellGeneDistances = dist.GetDistances(result)
#CellGeneDistances = dist.GetDistances(result, features=result.featuresCoordinates.index.tolist()[:10], cells=result.cellsCoordinates.index.tolist()[:50], dims=range(30))

print(CellGeneDistances)

# Optionally write distance results to CSV file
#CellGeneDistances.to_csv('Results_csv/CellGeneDistances.csv')
#CellGeneDistances.iloc[:15, :100].to_csv('Results_csv/test1.csv')












