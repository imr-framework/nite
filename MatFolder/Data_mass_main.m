%% This file describes the process of NITE DL recon prep 
% Author: Sairam Geethanath
% Input - Thermal simulations with tissue properties
% Output - Variables for training and testing ready for ELM
addpath(genpath('.'));
%% Get data ready in the right shape and range
Ns =4; %Number of temperature sensors at the surface
dT = 0.1; %C resolution required
datapath = 'MatFolder\Training_new.mat';

load(datapath);
Data = x_Data;
Temperatures =(squeeze(x_Data(:,end))); %0.01C round off

Np = size(Data,1)./Ns; %has to be an integer
[Data_scale,Tp_label,Tp_label_vec2] = Data_prep4ML(Data,Temperatures,Np,Ns,dT,'Rescale');

%% Data augment to enable DL - using rand noise for now
%% Data augment to enable DL - using rand noise for now
%compensate for lower number of examples for higher temperature values in
%training
counts = hist(Tp_label, length(Tp_label_vec2));
counts_max = max(counts);
counts_min = min(counts);

Data_scale_temp = Data_scale;
Tp_label_temp = Tp_label;
label_counter = [];

histogram(Tp_label,41);

for m = 1:size(Data_scale,1)
    if(~ismember(Tp_label(m),label_counter))
            label_counter = cat(2, label_counter, Tp_label(m));
            repvec = counts_max - counts(Tp_label(m) + 1);
            rep_mat = repmat(squeeze(Data_scale_temp(m,:)),[repvec 1]);
            rep_mat = rep_mat + 0.03.*rand(size(rep_mat));
            rep_mat(rep_mat > 1) = 1;
            Tp_label_curr = repmat(Tp_label(m),[repvec 1]);
            Data_scale = cat(1, Data_scale, rep_mat);
            Tp_label = cat(1,Tp_label, Tp_label_curr);
            disp(length(label_counter));
    end
end
figure;
hold on;

histogram(Tp_label,41);




if(length(Tp_label) < 10e3) %required for an accuracy of 0.96
            augfact = 12;
             max_noise = 0.05; %std of noise allowed as compared to signal of 0 to 1
            req_samp =length(Tp_label).*augfact; %number of samples required 
            [x2,t2] = data_augment(Data_scale, Tp_label, req_samp, max_noise);
else
       x2 = Data_scale;
       t2 = Tp_label;
end
 
 %% Partition to training/test and save
 train_percent = 95;% %of data set for training
 %test =100-train;
 partition_save(x2,t2,train_percent);
 

 
 %% Perform training and testing on Python
 %***********************************************************************
 %*********************************************************************** 

 %% Obtain the mat files from testing - labels and convert them back to temperature values
 