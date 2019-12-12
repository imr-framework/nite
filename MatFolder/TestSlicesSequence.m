%% Run this script to have the data prepared for all the brain slices (Train and Test): x_Data (HUGO) and all the t_Data (HUGO, NELLY)

cd MatFolder\

slice_TRAIN_HUGO_z70; %training slice (from HUGO)
fprintf('slice_TRAIN_HUGO_z70 DONE \n')
cd MatFolder\

slice_TEST_HUGO_z60; %test slice HUGO 1
fprintf('slice_TEST_HUGO_z60 DONE \n')
cd MatFolder\

slice_TEST_HUGO_z50; %test slice HUGO 2
fprintf('slice_TEST_HUGO_z50 DONE \n')
cd MatFolder\

Test_Slice_NELLY; %test slice NELLY 1
fprintf('Test_Slice_NELLY DONE \n')
cd MatFolder\

slice_TEST_HUGO_z90; %test slice HUGO 3
fprintf('slice_TEST_HUGO_z90 DONE \n')
cd MatFolder\

slice_TEST_HUGO_z80; %test slice HUGO 4
fprintf('slice_TEST_HUGO_z80 DONE \n')
cd MatFolder\

Test_Slice_NELLY_newSlice_z70; %test slice NELLY 2
fprintf('Test_Slice_NELLY_newSlice_z70 DONE \n')
cd MatFolder\

Test_Slice_NELLY_newSlice_z60; %test slice NELLY 3
fprintf('Test_Slice_NELLY_newSlice_z60 DONE \n')
cd MatFolder\

Test_Slice_NELLY_newSlice; %test slice NELLY 4
fprintf('Test_Slice_NELLY_newSlice DONE \n')
cd MatFolder\
