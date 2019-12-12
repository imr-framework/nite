function [xa,ta] = data_augment(x, t, req_samp, max_noise)
%% This needs better coding - more interpolation based rather than noise based
%% Create variables
xa = zeros(req_samp, size(x,2));
ta = zeros(req_samp, size(t,2));
ind = randperm(req_samp);
samp_fact = req_samp./size(x,1); %has to be an integer



%% Generate random variation of each sample
for k = 1:size(x,1)
    samp = squeeze(x(k,:)); tmp = repmat(samp, [samp_fact-1, 1]);
    label = squeeze(t(k,:));
    aug_op = rand(length(samp), samp_fact-1);
    samp_a = (rand(1).*max_noise.*aug_op).' +tmp;
    samp_a = cat(1, samp,  samp_a);
    ind_get = ((k-1)*samp_fact+1):(k*samp_fact);
    ind_use = ind(ind_get);
    xa(ind_use,:) = samp_a;
    ta(ind_use) = label;
end
