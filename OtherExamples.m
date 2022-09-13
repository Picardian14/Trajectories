% THIS FILE INCLUDES THE TIME CHANGING COMPARISON AND SWITCH OF DIMENSION
% SPACE EXAMPLES. SHOULD BE RUNNED AFTER EXECUTING N3DATAEXAMPLE.M


%% EPOCHED IN FRAGMENTS OF 100 AND 50
% I put pile them in the last dimension as if they were new trajectories
tsdata_epoched = cat(3, tsdata_success(:, 1:100, :),tsdata_success(:, 101:200, :));
data_reshape = tsdata_epoched;
n_rep = 512;
n_clust = 16;
start_clust = 2;
clusts = {};
clusts{end+1} = {};
for k=start_clust:n_clust
    [id, clust, distortion, silhouette, connectivity] = traj_kmeans(data_reshape, k, n_rep);
    
    for i = 1:k
        mat_out = data_reshape - repmat(clust(:, :, i),[1 1 size(data_reshape,3)]);
        dist_out(:, i) = squeeze(sum(sum(sum(abs(mat_out)))));
    end
    D(1, k) = k;
    D(2, k) = distortion;%sum(dist_out);
    clusts{end+1} = {clust, silhouette};
    C = mean(connectivity,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
    tmp = 0;
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            tmp = tmp + (4 * (C(i,j) - 1/2)^2);
        end
    end
    rho(k) = 1/((size(C,2))^2) * tmp;
end
figure
hold on;
plot(D(1,start_clust:end),D(2,start_clust:end));    
plot(D(1,start_clust:end),D(2,start_clust:end), 'or');
hold off
title("Elbow")
figure
hold on;
ks = start_clust:n_clust;
plot(ks,rho(start_clust:end));    
plot(ks,rho(start_clust:end), 'or');
hold off
title("Dispersion")
%%


k=2;
n_rep = 2048;
%[id_e, clust_e, distortion_e, silhouette_e, connectivity_e] = traj_kmeans(tsdata_epoched, k, n_rep);
C = mean(connectivity_e,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
tmp = 0;
for i = 1:size(C,1)
    for j = 1:size(C,2)
        tmp = tmp + (4 * (C(i,j) - 1/2)^2);
    end
end
rho = 1/((size(C,2))^2) * tmp;
figure
tiledlayout(1,2);    
colors = ["red", "black", "green", "blue"];
for roi = 1:size(clust,1)
    nexttile;
    hold on
    %for i_clust = 1:k
    %    plot(clust(roi, :, i_clust),'color', colors(1,i_clust));                        
    %end
    plot(clust_e(roi, :,1),'color', colors(1,1));                        
    plot(clust_e(roi, :,2),'color', colors(1,2));                        
                           
    title('Cluster '+string(roi))
    hold off     
end


classifications = zeros(k,2);
for i_clust=1:k
    classifications(i_clust, 1) = sum(id_e(1:NSUB)==i_clust)/NSUB;
    classifications(i_clust, 2) = sum(id_e(NSUB+1:NSUB*2)==i_clust)/NSUB;
end
classifications
%% Epoch of 50
%tsdata_epoched = cat(3, tsdata_pca(:, 1:50, :),tsdata_pca(:, 51:100, :),tsdata_pca(:, 101:150, :),tsdata_pca(:, 151:200, :));
%data_reshape = tsdata_epoched;
n_rep = 512;
n_clust = 16;
start_clust = 2;
clusts = {};
clusts{end+1} = {};
for k=start_clust:n_clust
    [id_e50, clust_e50, distortion_e50, silhouette_e50, connectivity_e50] = traj_kmeans(data_reshape, k, n_rep);
    
    for i = 1:k
        mat_out = data_reshape - repmat(clust_e50(:, :, i),[1 1 size(data_reshape,3)]);
        dist_out(:, i) = squeeze(sum(sum(sum(abs(mat_out)))));
    end
    D(1, k) = k;
    D(2, k) = distortion_e50;%sum(dist_out);
    clusts{end+1} = {clust_e50, silhouette_e50};
    C = mean(connectivity_e50,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
    tmp = 0;
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            tmp = tmp + (4 * (C(i,j) - 1/2)^2);
        end
    end
    rho(k) = 1/((size(C,2))^2) * tmp;
end
figure
hold on;
plot(D(1,start_clust:end),D(2,start_clust:end));    
plot(D(1,start_clust:end),D(2,start_clust:end), 'or');
hold off
title("Elbow")
figure
hold on;
ks = start_clust:n_clust;
plot(ks,rho(start_clust:end));    
plot(ks,rho(start_clust:end), 'or');
hold off
title("Dispersion")

%%

k=2;
n_rep = 2048;
%[id_e50, clust_e50, distortion_e50, silhouette_e50, connectivity_e50] = traj_kmeans(data_reshape, k, n_rep);
C = mean(connectivity_e50,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
tmp = 0;
for i = 1:size(C,1)
    for j = 1:size(C,2)
        tmp = tmp + (4 * (C(i,j) - 1/2)^2);
    end
end
rho = 1/((size(C,2))^2) * tmp;

figure
tiledlayout(1,2);    
colors = ["red", "black", "green", "blue"];
for roi = 1:size(clust,1)
    nexttile;
    hold on
    %for i_clust = 1:k
    %    plot(clust(roi, :, i_clust),'color', colors(1,i_clust));                        
    %end
    plot(clust_e50(roi, :,1),'color', colors(1,1));                        
    plot(clust_e50(roi, :,2),'color', colors(1,2));                        
                           
    title('Cluster '+string(roi))
    hold off     
end

classifications = zeros(k,2);
for i_clust=1:k
    classifications(i_clust, 1) = sum(id_e50(1:NSUB)==i_clust)/NSUB;
    classifications(i_clust, 2) = sum(id_e50(NSUB+1:NSUB*2)==i_clust)/NSUB;
end
classifications
%%
%%

% We have tN3o conditions, and we cant to see which ROIs characterize each
% one. We average accross subjects

mdata(1, :, :) = mean(Wdata, 3); % dimensionsxchannelsxtime
mdata(2, :, :) = mean(N3data, 3);
meandata = permute(mdata, [1 3 2]);

%% bad example: Here i Thiink its clustering two behaviours for each subject ??

%data = permute(tsdata, [3,2,1]); % Dimensions % time % channels
data = tsdata_pca_roi;
% Shape is 30x200x90
k = 2;
n_rep = 512;
[id, clust, distortion, silhouette, connectivity] = traj_kmeans(data, k, n_rep);
% Shape os clust is 30x200x4
% Clusters are levels of bold signal in 30 subjects
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
    %for i_clust = 1:k
    %    plot(clust(roi, :, i_clust),'color', colors(1,i_clust));                        
    %end
    plot(clust(roi, :,1),'color', colors(1,1));                        
    plot(clust(roi, :,2),'color', colors(1,2));                        
                           
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
    classifications(i_clust, 1) = sum(id(1:15)==i_clust)/NSUB;
    classifications(i_clust, 2) = sum(id(16:30)==i_clust)/NSUB;
end
classifications


%% Subject averaged clustering
meandata = permute(mdata ,[1 3 2]);
% Shape is 2 200 90
n_clust = 16;

n_rep = 6;
[id, clust] = traj_kmeans(meandata, n_clust, n_rep);
% Cluster is 2x200x16
%% Plotting

tiledlayout(4,4);   
colors = ["red", "black", "green", "blue"];
for i_clust= 1:n_clust   
    nexttile;    
    hold on
    plot(clust(2, :, i_clust),'color', "black", 'LineStyle','--');                        
    plot(clust(1, :, i_clust),'color', "red");                            
    title('Clust '+string(i_clust));
    hold off         
    
end

%{

for k=1:90
    msleep = squeeze(mean(NFCdataF, 1));
    mwake = squeeze(mean(WFCdataF, 1));
    sleep_r1 = msleep(:,k);
    wake_r1 = mwake(:,k);
    [h, p] = ttest(sleep_r1-wake_r1);
    if p>0.05
        disp(k)
    end
end


%}


