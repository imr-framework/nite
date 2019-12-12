# nite
Non Invasive Temperature Estimation for RF safety

Input data for these steps is available on Google Drive in CST_Data Folder.
If any of the links are broken, please open an issue.

1. [MATLAB] Data shaping after CST EM/Thermal simulation 

File: nite/MatFolder/TestSlicesSequence.m
Input: Google Drive links: CST_Data folder
Output: Mat files for each slices (training and test). Move output files to nite/MatFolder

2. [MATLAB] Data preparation for deep learning training 

File: nite/MatFolder/Data_mass_main.m 
Input: x_Data mat file 
Output: xtrain-4.mat, ytrain-4.mat, xtest-4.mat, ytest-4.mat 
Move output files to nite/PyFolder

3. [Python] Building and running the neural network for training 

File : nite/PyFolder/Train_Test.ipynb
Input: xtrain-4.mat, ytrain-4.mat, xtest-4.mat, ytest-4.mat  
Output: trained model (trained_model.h5)

4. [MATLAB] Data preparation for deep learning testing

File: nite/MatFolder/Data_mass_main_prep4test.m 
Input: t_Data mat files [changes with the considered slice]
Output: NITE_x_test_sl_HUGO4.mat with the variable x_test, NITE_x_test_sl_NELLY4.mat with the variable x_test [changes with the considered slice]
Move output files to nite/PyFolder

5. [Python] Testing slices with trained model  

File: nite/PyFolder/Training_Testing.ipynb
Input: NITE_x_test_sl_HUGO4.mat, NITE_x_test_sl_NELLY4.mat [changes with the considered slice]
Output: Reconstructed temperature maps (NITE_recon_HUGO_4-new.mat, NITE_recon_NELLY_4-new.mat)
Move output files to nite/MatFolder/


6. [MATLAB] Temperature maps plot and comparison between NITE map and CST map

Run nite/MatFolder/Data_mass_main.m again.
File: nite/MatFolder/Data_mass_main_testresvis.m
Input: NITE_x_test_sl_HUGO4.mat, NITE_x_test_sl_NELLY4.mat [changes with the considered slice] 
Output: Plot and comparison of temperature maps
