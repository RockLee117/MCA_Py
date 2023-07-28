Testing

How tests were performed:
-$$$ This MCA implementation has only been tested on Baron and Seger pancreas single-cell RNA-seq data sets provided in Baron et al. 2016 and Segerstolpe et al. 2016.
    https://storage.googleapis.com/cellid-cbl/BaronMatrix.rds
    https://storage.googleapis.com/cellid-cbl/BaronMetaData.rds
    https://storage.googleapis.com/cellid-cbl/SegerstolpeMatrix.rds
    https://storage.googleapis.com/cellid-cbl/SegerstolpeMetaData2.rds
-Using an RStudio project that uses the github link https://github.com/RausellLab/CelliD (CellID package), I added source() and sourceCpp() function calls in the R files to be able to manually call nested functions within RunMCA(), and call the RunMCA() located in my mca.R file with print statements to see results, the added source() calls, and write.csv statements to debug python implementation of RunMCA and to get the distance results in R. Look at the RFiles folder on how to setup the RStudio environment/files I used for testing.
-$$$$ Results using CellID in RStudio were considered to be correct then results from my python implementation were considered to be predictions
-(Goal was to mimic results from the R CellID package being performed when running MCA on same input data and distance results between the coordinates of the features and cells)
-Wrote python results of the feature & cell coordinates and the distance matrix to csv files
-Observed results in microsoft excel and values between R and python were relatively the same (may be difference by the 4th decimal). Note for both R and Python implementations, the resulting coordinates sign may change every time you run RunMCA but the absoulte values are the same. This results in the orientation of the coorindates to be different every time you run RunMCA but the distance between gene and cells are the same no matter the orientation
-To validate entire output, I performed mean absolute error on data considering the R output to be correct and the python output to be the predictions. Also, tested to see if the shapes (num rows and num columns) and the row and column names were the same. This is done in test.py. Follow comments in test.py to see how testing was done. Results are displayed in the .png files


How you can test yourself:
-In lines 17-22 in test.py, must put path to your csv outputs. test.py is setup assuming you write the coordinate outputs without the row names and the distance output to include the rownames
-To generate R output csv files, look at SettingUpRStudioForTestingCellId.txt to setup my R environment then in mca.R uncomment code that writes the coordinates to csv files and in RunMCATesting.R run the write.csv line after getting distance results.
-move the csv files to where your python code is located
-Before running RunMCA in python, make sure to follow the main ReadMe documentation on getting either the Baron or Seger datasets in readable format to be used in python
-To generate Python output csv files, copy dist.py, mca_steps.py, mca.py, and main.py into a folder with your input .csv files of the assay and row names. In main.py, point lines 25 and 27 to the input csv files. Then run main.py and uncomment the lines that write the coordinates and distances to csv files
-After adding paths in test.py, run test.py and it will output resuling mean absolute error values for the feature & gene coordinates and the distance results