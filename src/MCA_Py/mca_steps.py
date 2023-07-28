import numpy as np

"""
Computes fuzzy matrix
input: 
    - X: a 2D numpy array of the input dataset passed into RunMCA() and after preprocessing has been done
returns a list:
    - Z: numpy 2D array   
    - Dc: numpy vector (array)
"""
def MCAStep1(X):
    AM = X
    rmin = np.min(AM, axis=1) #get min val in each row
    rmax = np.max(AM, axis=1) #get max val in each row
    range_list = rmax - rmin 
    # subtract the min value for a row from each scalar in a row
    for index in range(len(rmin)):
        AM[index] -= rmin[index]
    # divide each value in a row vector by the row's range difference of (max - min) in a row
    for index in range(len(range_list)):
        AM[index] /= range_list[index]
        AM[index] = np.nan_to_num(AM[index]) # turns nan values to 0 in case where a number is divided by 0
    # fuzzy matrix stored in FM
    FM = np.concatenate((AM, 1 - AM), axis=0).astype('float64')
    AM = None
    total = np.sum(FM)
    colsum = np.sum(FM, axis=0)
    rowsum = np.sum(FM, axis=1)
    for index in range(len(rowsum)): # iterate the rows
        FM[index] /= np.sqrt(colsum)
        FM[index] = np.nan_to_num(FM[index]) # turns nan values to 0 in case where a number is divided by 0
    for index in range(len(rowsum)):
        FM[index] /= np.sqrt(rowsum[index])
        FM[index] = np.nan_to_num(FM[index]) 
    Dc = 1/(np.sqrt(rowsum/total))
    Dc = np.nan_to_num(Dc) 
    Dc = np.reshape(Dc, (np.shape(Dc)[0], 1))# turn Dc into a column vector 
    return {"Z": FM , 
            "Dc": Dc} 

"""
Computes coodinates
input:
    - Z: numpy 2D array from MCAStep1
    - V: numpy 2D array of Right singuar vectors from SVD computation
    - Dc: numpy vector (array) from MCAStep1
returns a list:
    - cellsCoordinates: numpy 2D array of cells coordinates where rows are cells and there are nmcs number of columns where columns represent dimensions
    - features: numpy 2D array of features coordinates where rows are features (genes) and there are nmcs number of columns where columns represent dimensions
"""
def MCAStep2(Z, V, Dc):
    AV = V 
    AZ = Z
    ADc = Dc
    Dc.reshape(1, np.shape(Dc)[0]) # turn column vector into a row vector
    FeaturesCoordinates = np.dot(AZ , AV)
    AZcol = np.shape(AZ)[1] # num columns in AZ
    AZ = None
    for index in range(np.shape(FeaturesCoordinates)[0]):
        FeaturesCoordinates[index] *= ADc[index]
    ADc = None
    cellsCoordinates = np.sqrt(AZcol) * AV
    return {"cellsCoordinates": cellsCoordinates,
            "featuresCoordinates": FeaturesCoordinates[ :int(np.shape(FeaturesCoordinates)[0]/2)] }


######### Personal Testing ########
# X = np.arange(1,26).reshape(5,5).astype('float64')
# result = MCAStep1(X) 
# print(result["Z"])
# print(result["Dc"])
 
# Z = np.arange(1,26).reshape(5,5).astype('float64')
# V = np.arange(1,26).reshape(5,5).astype('float64')
# Dc = np.arange(1,6).astype('float64').reshape(5,1)

# X = MCAStep2(Z,V,Dc)
# print(X["cellsCoordinates"])
# print(X["featuresCoordinates"])