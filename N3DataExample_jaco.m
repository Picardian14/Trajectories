
%Read the empirical data 

clear all;
close all;

load('DataSleepW_N3.mat'); % your data here, you need ts + structural connectivity
% se normaliza la conectividad estructural para que no sea mayor a 0.2 con el fin de que el acople no sea muy grande
C=SC/max(max(SC))*0.2;

N=90;
NPARCELLS = N;
NSUB=15;
Tmax=200;
Tepoch=20;
Toverlap=10; % a implementar despues
indexsub=1:NSUB;
indexepoch = 1:((Tmax-Tepoch)/Toverlap+1);
Isubdiag = find(tril(ones(N),-1));

params = [];
params.flp = 0.04; % Could try 0.01 0.1
params.fhi = 0.07;
params.TR = 2;


% Separing, N3, W, and N3+W
for nsub=indexsub
        
    TS = TS_W{1, nsub}/norm(TS_W{1, nsub});    
    tsdata=TS(:,1:Tmax) ; 
    tsdataF = permute(filter_bold(tsdata', params.flp, params.fhi, params.TR), [2 1 3]);
    tsdataF = tsdataF/norm(tsdataF); 
    
    for nepoch=indexepoch
        data_epochW(:,:,nepoch,nsub)=tsdataF(:,(1+(nepoch-1)*Toverlap):((nepoch-1)*Toverlap+Tepoch));
    end


    TS = TS_N3{1, nsub}/norm(TS_N3{1, nsub});    
    tsdata=TS(:,1:Tmax) ; 
    tsdataF = permute(filter_bold(tsdata', params.flp, params.fhi, params.TR), [2 1 3]);
    tsdataF = tsdataF/norm(tsdataF); 
    
    for nepoch=indexepoch
        data_epochN3(:,:,nepoch,nsub)=tsdataF(:,(1+(nepoch-1)*Toverlap):((nepoch-1)*Toverlap+Tepoch));
    end

end


%% Zscore 

% concatenate all data to compute norm 

all_data = cat(3,data_epochN3(:,:,:),data_epochW(:,:,:));

mean_value = mean(all_data(:,:),2);
std_value  = std(all_data(:,:),[],2);

all_data = (all_data-repmat(mean_value,[1,size(all_data,2),size(all_data,3)]))./repmat(std_value,[1,size(all_data,2),size(all_data,3)]);


%% baseline

all_data = all_data - repmat(all_data(:,1,:),[1 size(all_data,2) 1]);

 %% absolute value
 
 all_data = abs(all_data);

%%

%  tsdata = N3+W
for nsub=indexsub
    tsdata(:,:,nsub+NSUB)=TS_N3{1,nsub}(:,1:Tmax);    
    tsdataF(:,:,nsub+NSUB) = permute(filter_bold(tsdata(:, :,nsub+NSUB)', params.flp, params.fhi, params.TR), [2 1 3]);    
    tsdataF(:, :, nsub+NSUB) = tsdataF(:, :, nsub+NSUB)/norm(tsdataF(:, :, nsub));
end
%% Remove 2 for validation
val(:, :, 1) = tsdataF(:, :, 1);
val(:, :, 2) = tsdataF(:, :, NSUB+1);
tsdataF(:, :, 1) = [];
tsdataF(:, :, NSUB+1) = [];
NSUB = 14;
indexsub=1:NSUB;
%%
success_rois = SearchSignificantROI(TS_W, TS_N3, indexsub);
% When i have all subjects i find 88 and 90 as significantly different
% regions
% success_rois = [88 90];
%% good example: Here i should be clustering Whole brains of subjects into clusters
% Each trajectory would be a brain region dimensional time lapse of each
% subject whose condition can be wake or sleep
%tsdata_success = tsdataF(success_rois, :, :);
data_reshape = all_data;
k=4;
n_rep = 2048;
[id, clust, distortion, silhouette, connectivity, clust_size] = traj_kmeans(data_reshape, k, n_rep);
C = mean(connectivity,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
tmp = 0;
for i = 1:size(C,1)
    for j = 1:size(C,2)
        tmp = tmp + (4 * (C(i,j) - 1/2)^2);
    end
end
rho = 1/((size(C,2))^2) * tmp;

figure
tiledlayout(5,3);    
colors = ["red", "black", "green", "blue"];
for roi = 16:30
    nexttile;
    hold on
    for i_clust = 1:k
        plot(clust(roi, :, i_clust),'color', colors(1,i_clust));                        
    end
%     plot(clust(roi, :,1),'color', colors(1,1));                        
%     plot(clust(roi, :,2),'color', colors(1,2));                        
%     plot(clust(roi, :,3),'color', colors(1,3));   
%     plot(clust(roi, :,4),'color', colors(1,4));   
    title('ROI PCA'+string(roi))
    hold off     
end
figure
tiledlayout(1, size(clust_size, 1));
for t=1:size(clust_size, 1)    
    nexttile;
    histogram(clust_size(t, :));   
end
std(clust_size');

classifications = zeros(k,2);
for i_clust=1:k
    classifications(i_clust, 1) = sum(id(1:150)==i_clust)/(NSUB*10);
    classifications(i_clust, 2) = sum(id(151:300)==i_clust)/(NSUB*10);
end
classifications
%%
%val = val(success_rois, :, :)
for i=1:k
    mat = val - repmat(clust(:,:,i),[1 1 size(val,3)]);
    dist(:,i) = squeeze(sum(sum(abs(mat))));
end 
[a_val,id_val] = min(dist'); 


