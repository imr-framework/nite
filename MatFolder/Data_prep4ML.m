function [Data_scale,Tp_label,Tp_label_vec] = Data_prep4ML(Data,Temperatures,Np,Ns,dT,feature_scale)
%Create dielectric properties matrix 
% Np = 115; %m - number of points/examples ~1e4
% ndist =11;
% Ns = 4;
%% Feature scaling - Obtain measures for scaling
Dmax = max(Data,[],1); Dmax_table = repmat(Dmax, [size(Data,1) 1]);
Dmin = min(Data,[],1);  Dmin_table = repmat(Dmin, [size(Data,1) 1]);
Davg = mean(Data,1);  Davg_table = repmat(Davg, [size(Data,1) 1]);
Dstd = std(Data,[],1);  Dstd_table = repmat(Dstd, [size(Data,1) 1]);
% feature_scale = 'Rescale';
disp(['Chosen feature scaling method: ', feature_scale]);
switch feature_scale
    case 'None'
    Data_scale = Data;
    case 'Rescale'
    Data_scale = (Data - Dmin_table)./(Dmax_table - Dmin_table + 0.06); %Dmax and Dmin are the same causing redundancy here - remove the 0.1 in the next iteration
    case 'Normalize'
    Data_scale = (Data - Davg_table)./(Dmax_table - Dmin_table + 0.06);    
    case 'Standardize'
    Data_scale = (Data - Davg_table)./(Dstd_table);    
     %case 'Unit_norm'
     %Data_scale = (Data)./(Dnorm);        
end
    

%% Prepare Data_scale 
if(length(Temperatures) > 1) %training case
    Data_scale = squeeze(Data_scale(:,1:end-1));
end

 Data_scale = reshape(Data_scale, [Np Ns*35]); %35 features 
%% Handle y to convert regression to classification
Tp = Temperatures(1:Ns:end);
Tp_label_vec = min(Tp):dT:max(Tp);  %Looking for a 0.1C resolution 
Tp_label = zeros(size(Tp));
for p = 1:length(Tp)
    dist_vec = abs(Tp_label_vec - Tp(p));
    [~, label] = min(dist_vec);
    Tp_label(p) = label - 1; %to be consistent with Python convention of 0 to N-1
end

