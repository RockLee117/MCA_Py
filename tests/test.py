import numpy as np
import pandas as pd
# File is responsible for reading the outputs of the coordinates and distances from my python implementation and R's CellId
# Then calculating mean absolute error for coordinates and distances to test accuracy of the model

# mean absolute error
def MAE(true, predictions):
    true, predictions = np.array(true), np.array(predictions)
    return np.mean(np.abs(true - predictions))

# NOTE: A MAE of 0 indicates that the model's predictions are perfectly accurate, and there are no errors between the predicted and actual values. Closer MAE value is to 0 the better

# Testing MCA Results using Baron dataset
def test():
    
    # Baron dataset results
    R_cells_coordinates = pd.read_csv('ROutputData/R_cellsCoordinates.csv')
    R_features_coordinates = pd.read_csv('ROutputData/R_featuresCoordinates.csv')
    R_distances = pd.read_csv('ROutputData/R_FullBaron_DISTANCES.csv')
    Py_cells_coordinates = pd.read_csv('../Results_csv/cellsCoordinates.csv')
    Py_features_coordinates = pd.read_csv('../Results_csv/featuresCoordinates.csv')
    Py_distances = pd.read_csv('../Results_csv/CellGeneDistances.csv')
    
    # Seger dataset results
    # R_cells_coordinates = pd.read_csv('ROutputData/Seger_R_cellsCoordinates.csv')
    # R_features_coordinates = pd.read_csv('ROutputData/Seger_R_featuresCoordinates.csv')
    # R_distances = pd.read_csv('ROutputData/R_Seger_distances.csv')
    # Py_cells_coordinates = pd.read_csv('../Results_csv/Seger_cellsCoordinates.csv')
    # Py_features_coordinates = pd.read_csv('../Results_csv/Seger_featuresCoordinates.csv')
    # Py_distances = pd.read_csv('../Results_csv/Seger_CellGeneDistances.csv')
    
    # exclude first column b/c those are the row names (the data frames containing distance values were written to csv files INCLUDING THE ROW NAMES)
    rownames = R_distances.iloc[:, 0].tolist()
    R_distances.index = rownames
    R_distances = R_distances.iloc[:, 1:]
    
    rownames = Py_distances.iloc[:, 0].tolist()
    Py_distances.index = rownames
    Py_distances = Py_distances.iloc[:, 1:]
    
    # For coordinates, go column by column and perform the mean absolute error calculation
    # Then average the mean absolute error generated for each column
    
    # cell coordinates test
    if R_cells_coordinates.shape != Py_cells_coordinates.shape:
        print('Error: R_cells_coordinates and Py_cells_coordinates do not have the same shape')
        return
    print('Cell coordinates are same shape')

    sum_MAE = 0.0 # will hold sum of the mean absoulte error for each column
    num_cols = R_cells_coordinates.shape[1] 
    
    for i in range(R_cells_coordinates.shape[1]):
        # NOTE: must take absolute value before performing MAE because of varying output sign (positive or negative)  for the coordinates everytime you run RunMCA(). (but value ends up being same)
        mae_res = MAE(np.abs(R_cells_coordinates.iloc[:, i]), np.abs(Py_cells_coordinates.iloc[:, i]))
        sum_MAE += mae_res
    print('Average of the mean absolute error for CELL coordinates is: ', sum_MAE/num_cols)
        
        
    # feature coordinate test
    if R_features_coordinates.shape != Py_features_coordinates.shape:
        print('Error: R_features_coordinates and Py_features_coordinates do not have the same shape')
        return
    print('Features coordinates are same shape')

    sum_MAE = 0.0 # will hold sum of the mean absoulte error for each column
    num_cols = R_features_coordinates.shape[1] 
    
    for i in range(R_features_coordinates.shape[1]):
        # cell coordinates test
        # NOTE: must take absolute value before performing MAE because of varying output sign (positive or negative)  for the coordinates everytime you run RunMCA(). (but value ends up being same)
        mae_res = MAE(np.abs(R_features_coordinates.iloc[:, i]), np.abs(Py_features_coordinates.iloc[:, i]))
        sum_MAE += mae_res
    print('Average of the mean absolute error for FEATURES coordinates is: ', sum_MAE/num_cols)
    
    # distance calculation MAE
    if R_distances.shape != Py_distances.shape:
        print('Error: distance calculations are not the same shape')
        return
    print('Distance calculations are the same shape')
    
    if set(R_distances.index) != set(Py_distances.index):
        print('Error: R_distances and Py_distances do not have the same row names')
        return
    print('R_distances and Py_distances have the same row names')
    
    if set(R_distances.columns) != set(Py_distances.columns):
        print('Error: R_distances and Py_distances do not have the same column names')
        return
    print('R_distances and Py_distances have the same column names')
    
    sum_MAE = 0.0 # will hold sum of the mean absoulte error for each column
    num_rows = R_distances.shape[0] 
    # average the MAE of all the rows b/c we want to measure the accuracy of distance of Gene X from all of the cells (columns)
    for i in range(R_distances.shape[0]):
        mae_res = MAE(R_distances.iloc[i, :], Py_distances.iloc[i, :])
        sum_MAE += mae_res
    print('Average of the mean absolute error for DISTANCES is: ', sum_MAE/num_rows)
        
test()