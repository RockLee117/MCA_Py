<h1 align="center">MCA_Py</h1>

## Overview

Based on the MCA performed in the Cellid R package: https://github.com/RausellLab/CelliD

MCA_Py performs multiple correspondence analysis, in Python, on input 2D NUMERICAL, Tabular Data. By calling the function RunMCA, the output are 2 pandas data frames, one being the coordinates of input dimension x (rows from the input matrix) and the other data frame being coordinates for input dimension y (columns from the input matrix) in a j dimensional space. By default, the output coordinates contain 50 dimensions where the columns in the coordinate data frames are the ith dimension. <br>
In addition, you are able to find the distance between the two coordinate data frames by calling the function GetDistances. The output will be a pandas data frame containing values representing the distance of row x to column y based off of the coordinates.

IMPORTANT: This MCA implementation has ONLY been tested on Baron and Seger pancreas single-cell RNA-seq data sets provided in <a href="https://www.sciencedirect.com/science/article/pii/S2405471216302666?via%3Dihub">Baron et al. 2016</a> and <a href="https://www.sciencedirect.com/science/article/pii/S1550413116304363?via%3Dihub">Segerstolpe et al. 2016</a> where genes are rows, columns are cells, and each value represents the gene expression of gene x in cell y. <br>
https://storage.googleapis.com/cellid-cbl/BaronMatrix.rds <br>
https://storage.googleapis.com/cellid-cbl/BaronMetaData.rds <br>
https://storage.googleapis.com/cellid-cbl/SegerstolpeMatrix.rds <br>
https://storage.googleapis.com/cellid-cbl/SegerstolpeMetaData2.rds <br>

## Installation
Assuming python is installed on your system along with pip, MCA can be installed running: <br>
pip install MCA_Py

## How to Use Step by Step

1. First, store your input data into a pandas data frame with the appropriate row and column names. If trying to use the Baron or Seger datasets, follow steps in loadAssay.txt in the doc folder. Note the input data is expected to be numerical values and not contain any strings.
2. Using the imported MCA package, call RunMCA and pass in your pandas data frame. The resulting coordinate data frames will be stored in an object along with a vector of the singular values resulting from the singular value decomposition performed during the MCA process. IMPORTANT NOTE: the sign of the values in the output coordinates may vary everytime you call RunMCA with the same input data. This causes the orientation of the coordinates to change, BUT the distance between the coordinates remains the same. For example, one time a value may be 4 and another time may be -4.
3. Optionally, you can graph the coordinates using matplotlib on a scatterplot, but must pick only 2 of the (default) 50 dimensions. Choose one column in both of the coordinate data frames to be the x and another column to be the y. Ideally plot the first column as the x and the second column as the y.
4. Optionally, you can find the distances between the coordinates of the 2 data frames by calling GetDistances and passing in the object returned from RunMCA. The result will be a pandas data frame containing the distances.

(Feel free to refer to main.py in the MCA_Py github to see an example of performing MCA and finding the distances between coordinates using the SOURCE code files, dist.py, mca_steps.py and mca.py)

## Example

Here is an example of using the MCA_Py package on the Baron pancreas single-cell RNA-seq data set: <br>

import pandas as pd <br>
import numpy as np <br>
import matplotlib.pyplot as plt <br>
import MCA_Py, csv

<br># Step 1: Store input into pandas dataframe<br>
<br># (Assuming using Baron dataset and looked at loadAssay.txt in the doc folder to see how to get necessary csv files)
<br># Assay of dataset 
<br>assay_df = pd.read_csv('Baron_assay_data.csv')  # read_csv contains string path to csv file containing the values in the dataset and column names
<br>
<br># row names of assay 
<br>row_names = pd.read_csv('Baron_assay_rownames.csv') # csv contains row names for the datset
<br>
<br># set row names to data frame's row names 
<br>assay_df.index = np.ndarray.flatten(row_names.values)

<img src="https://github.com/RockLee117/Images/blob/main/MCAOutputImages/input.png" width=500>
<br>
<br>

<p># Step2: Perform MCA on input data stored in pandas data frame</p>

"""
<br>
result is an object containing these members: <br>
 - featuresCoordinates: pandas data frame where rows are genes and columns are nmcs (default value for nmcs is 50) different dimensions <br>
 - cellsCoordinates: pandas data frame where rows are cells and columns are nmcs (default value for nmcs is 50) different dimensions <br>
 - stdev: numpy array containing the singular values created during MCA when performing Singular Value Decomposition <br>
"""
<br>
result = MCA_Py.RunMCA(assay_df) <br>

#Optional parameters for RunMCA(): <br>
- nmcs: (number of output dimensions to compute in coordinate matrices). Default value is 50 <br>
- features: string list specifying which rows to use by listing the row names. Default value is None so ends up using all the rows if not specified

<img src="https://github.com/RockLee117/Images/blob/main/MCAOutputImages/mca_output.png" width=500> <br>
<img src="https://github.com/RockLee117/Images/blob/main/MCAOutputImages/mca_singularvals.png" width=500>
<br>

<br># Step3: Optionally plot coordinates. Here, the first 2 dimensions/columns are being plotted where the first column is the x and the second is the y in the scatterplot<br>
<br>plt.scatter(result.featuresCoordinates.iloc[:, 0], result.featuresCoordinates.iloc[:, 1], label='Genes(features) Coordinates', c='blue', marker='x')
<br>plt.scatter(result.cellsCoordinates.iloc[:, 0], result.cellsCoordinates.iloc[:, 1], label='Cells Coordinates', c='red', marker='o', alpha=0.5)
<br># Add labels and legend
<br>plt.xlabel('MCA_1')
<br>plt.ylabel('MCA_2')
<br>plt.title('MCA Scatter Plot')
<br>plt.legend()
<br># Show the plot
<br>plt.show() <br>

<img src="https://github.com/RockLee117/Images/blob/main/Python_EntireBaron.png" width=500>

<br># Step4: Optionally calculate the distances
<br>CellGeneDistances = MCA_Py.GetDistances(result)

#Optional parameters for GetDistances:<br>
 - dims: specifies a range of which dimensions/columns to use in the input coordinate data frames. Important: if a dims value was specified when calling RunMCA, then must specify a range value for dims when calling GetDistances<br>
 - features: String vector of feature names (rows) to subset from feature coordinates. If not specified will use all features 
- cells: String vector of cell names (rows) to subset from cell coordinates. If not specified will use all cells 
#

<img src="https://github.com/RockLee117/Images/blob/main/MCAOutputImages/dist.png" width=500>