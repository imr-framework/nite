%% This file describes the process of NITE DL recon prep 
% Author: Sairam Geethanath
% Input - labels from ELM test
% Output - Correlation with ground truth labels
addpath(genpath('.'));
%datapath = 'C:\Users\Julie\Desktop\NITE_May2019\Matlab_Code_2';
%matname = 'NITE_x_test_sl_NELLY.mat';

datapath0 = 'MatFolder\';
matname = 'Test_HUGOz90.mat'; % to be changed for each test slice 
%%
Ns = 4;
dT =0.1;
datapath = fullfile(datapath0,matname); 
load(datapath);
Data = t_Data_z90;%changes as per the data set - IMPORTANT !!!!
Temperatures = 0;%x_Data(:,end); to distinguish from training
Np = size(Data,1)./Ns; %has to be an integer
[Data_scale,Tp_label,Tp_label_vec2] = Data_prep4ML(Data,Temperatures,Np,Ns,dT,'Rescale');
x_test4 = Data_scale.';
cdir = pwd;
cd(datapath0);
save('NITE_x_test_sl_HUGO4','x_test4');
cd(cdir);