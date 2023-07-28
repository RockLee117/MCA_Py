<h1 align="center">MCA</h1>

## Overview

Based on the MCA performed in the Cellid package: https://github.com/RausellLab/CelliD

MCA performs multiple correspondence analysis in Python on input 2D NUMERICAL Tabular Data. Output are 2 pandas data frames, one being the coordinates of input dimension x (rows from the input matrix) and the other data frame being coordinates for input dimension y (columns from the input matrix) in a j dimensional space (by default, output coordinates contain 50 dimensions). The columns in the coordinate data frames are the ith dimension. <br>
In addition, you are able to find the distance between the two coordinate matrices by calling GetDistances. The output will be a pandas data frame containing values representing distance of row x to to column y (based off the coordinates).

IMPORTANT: This MCA implementation has ONLY been tested on Baron and Seger pancreas single-cell RNA-seq data sets provided in Baron et al. 2016 and Segerstolpe et al. 2016 where genes are rows and columns are cells. Each value represents the gene expression of gene x in cell y. <br>
    https://storage.googleapis.com/cellid-cbl/BaronMatrix.rds
    https://storage.googleapis.com/cellid-cbl/BaronMetaData.rds
    https://storage.googleapis.com/cellid-cbl/SegerstolpeMatrix.rds
    https://storage.googleapis.com/cellid-cbl/SegerstolpeMetaData2.rds

## Installation
Assuming python is installed on your system along with pip, MCA can be installed running: <br>
pip install sdvssdfsdgddfg

## How to Use Step by Step
-Demonstration using main.py,- Step by step order of calling functions
-show how to optionally write results to csv and graph coordinates in scatterplot

1. First, store your input data into a pandas data frame with the appropriate row and column names. If trying to use the Baron or Seger datasets, follow steps in loadAssay.txt
2. Using the imported MCA package, call RunMCA and pass in your pandas data frame. The resulting coordinate data frames will be stored in an object along with a vector of the singular values resulting from the singular value decomposition performed during the MCA process.
3. Optionally, you can find the distances between the coordinates of the 2 coordinate data frames by calling GetDistances and passing in the object returned from RunMCA. The result will be a pandas data frame.

Here is an example of using the MCA package: <br>





## Results
-Results using Baron Dataset only using the genes that fall in the HgProteinCodingGenes and Seger dataset and compare them to R output (visually show images of the values with row and column names)
