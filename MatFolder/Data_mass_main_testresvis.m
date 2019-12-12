addpath(genpath('.'));
datapath = 'MatlabFolder\';
%%
Ns = 4;
dT =0.1;
%% Post running the DL code
        load(fullfile(datapath,'NITE_recon_NELLY_2.mat')); %Loads estimated Y_hat - "labels" - change path according to the test slice to be evaluated after DL testing 
        labels = ynew2;
        %load('NITE_t_test'); %Loads GT labels - "t_test"
        %load('label2temp.mat'); %Loads label2temperature vector - "Tp_label_vec"
        M_test = M_z70_N;
        Tp_GT = squeeze(M_test(:,3));
        X = squeeze(M_test(:,1));
        Y = squeeze(M_test(:,2));
        t_test = zeros(size(Tp_GT));
        
        %Bin GT temperature to compare against bins
        for p = 1:length(t_test)
                dist_vec = abs(Tp_label_vec2 - Tp_GT(p));
                [~, label] = min(dist_vec);
                t_test(p) = label -1; %to be consistent with Python convention of 0 to N-1
        end
        
% Convert labels to temperature
            Tp_GT = label2temp(t_test, Tp_label_vec2);
            %labels=double(labels);
            Tp_NITE = label2temp(labels, Tp_label_vec2);
             
% Display
%             figure;
%             plot(t_test,labels,'r*');
%             xlabel('GT labels'); ylabel('Estimated labels');
%             grid on;
          %  Classification output
            figure; title ('Training temperature correlations');
            plot(Tp_GT,Tp_NITE.','r*');
            xlabel('CST (^0C)'); ylabel('NITE (^0C)');
            grid on;

          % Temperature voxels output
            figure;
            plot(Tp_GT,'k*')
            grid on; hold on;
            plot(Tp_NITE,'bo');legend('CST','NITE');
            xlabel('Sample #'); 
            ylabel('Temperature  (^0C)');
           
           % Spatial maps
           %load(fullfile(datapath, 'XY_coord')); 
           figure;
           scatter(X, Y, 20, Tp_NITE, 'filled');
     
            colorbar
            title('Temperature map after 60s: Test slice NITE [°C]')
            set(gca,'clim',[37 42.5])
            set(gcf,'Position',[600 0 550 600])


            %load(fullfile(datapath, 'XY_coord')); 
            figure;
            scatter(X, Y, 20, Tp_GT, 'filled');
           

            colorbar
            title('Temperature map after 60s: Test slice GT [°C]')
            set(gca,'clim',[37 42.5])
            set(gcf,'Position',[600 0 550 600])


%%%%%%%%%%%%%




