function  partition_save(x,t,train_percent)
%%
num_classes = 41;%fixed for NITE
samples_train = ceil(size(x,1).*train_percent./100); 

rand_vec = randperm(size(x,1));


xtrain2 = squeeze(x(rand_vec(1:samples_train),:)).';
xtest2  = squeeze(x(rand_vec(samples_train+1:end),:)).';
ytrain_temp2  = squeeze(t(rand_vec(1:samples_train),:)).'; %save this guy to test on training - classes 0 - 46
ytest_temp2 = squeeze(t(rand_vec(samples_train+1:end),:)).';%save this for comparing with the validation

%% create one hot vector for classification
ytrain2 = zeros(num_classes, size(ytrain_temp2,2));
for k = 1:size(ytrain2,2)
    m = ytrain_temp2(k) +1;
    ytrain2(m,k) = 1;
end

ytest2 = zeros(num_classes, size(ytest_temp2,2));
for k = 1:size(ytest2,2)
    m = ytest_temp2(k) +1;
    ytest2(m,k) = 1;
end





%% Save variables as mat files
% save('NITE_x_train', 'x_train');
% save('NITE_t_train', 't_train')
% save('NITE_x_test', 'x_test')
% save('NITE_t_test','t_test')

save('xtrain-4', 'xtrain2');
save('ytrain-4', 'ytrain2')
save('xtest-4', 'xtest2')
save('ytest-4', 'ytest2')
save('ytrain_store-1', 'ytrain_temp2')
save('ytest_store-1', 'ytest_temp2')
