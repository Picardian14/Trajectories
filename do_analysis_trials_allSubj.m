clear all
close all

%% Parameters %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ex_subj=[4 10 16 18 21]; % excluded subjects
Condition='Active'; %;'Passive'
n_comp = [1:2]; % PCA comp that are used for clustering
n_cluster = 2; % how many clusters
thresh_SNR=0; %cluster data with SNR >0
mid_intensity=3; %Intensity levels for checking response in clusters of equal intensity 

%% Prepare CSD matrix with all subj --> completed and saved as data_csd_allsubj.mat

data_csd = [];
intensity = []; 
response = [];
Sbj_id = [];

for sbj=1:25
    if ismember(sbj,ex_subj)==0 % continue if current subj is not excluded
        
        if sbj<10
          
            data_csd_aux = load (['C:\Users\tal.malkinson\Documents\Trajectory Kmeans\Data\data_CSD_Sub0' num2str(sbj) '_' Condition '.mat']);
            load (['C:\Users\tal.malkinson\Documents\Trajectory Kmeans\Data\data_ref_Sub0' num2str(sbj) '_' Condition '.mat']);
        else
            data_csd_aux = load (['data_CSD_Sub' num2str(sbj) '_' Condition '.mat']);
            load (['data_ref_Sub' num2str(sbj) '_' Condition '.mat']);
        end
        
        data_csd = cat(3,data_csd , data_csd_aux.data_csd);
        intensity = [intensity ; data_ref.trialinfo(:,2)];
        response = [response ; data_ref.trialinfo(:,6)];
        Sbj_id = [Sbj_id ; ones(size(data_csd_aux.data_csd,3),1).* sbj];
        
        clearvars data_csd_aux data_ref
    end
end

%% 
% cd C:\Users\tal.malkinson\Documents\Trajectory Kmeans\Data
% 
% load C:\Users\tal.malkinson\Documents\Trajectory Kmeans\Data\data_csd_allsubj.mat
% 
% cd C:\Users\tal.malkinson\Documents\Trajectory Kmeans

% zscore data

data_csd = zscore(data_csd,[],2);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cluster trials  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%




%%% keep SNR above 0 

trials = intensity>thresh_SNR;

data = data_csd(:,:,trials);


%%%% pca decomposition

disp(' ******* Computing PCA *******')

[coeff,score,latent,tsq, explained, mu] = pca(mean(data,3)');


figure
imagesc(mean(data,3))
title('original data')

for i = 1:size(data,3)
   
%    data_pca(:,:,i) = data(:,:,i)'*score(:,1:10);
    data_pca(:,:,i) = coeff(:,n_comp)'*data(:,:,i);
    
end


%data_pca = permute(data_pca,[2 1 3]);


figure
imagesc(mean(data_pca,3))
title('data in PCA space')


%%% reconstruct data

for i = 1:size(data_pca,3)
   
    data_reconstructed(:,:,i) = coeff(:,n_comp)*data_pca(:,:,i);
    
end

figure
imagesc(mean(data_reconstructed,3))
title('data reconstructed from PCA space')

figure
ecdf(explained)
title('cumulative explained variance by PCA components')



%%% behavior
response = response(trials);
intensity = intensity(trials);

%% Clustering


disp(' ******** CLUSTERING *******')
[id,cluster_out] = traj_kmeans(data_pca, n_cluster, 50);

figure('Name','All Intensities - GFP of clusters on PCA space','NumberTitle','off');


subplot(221)
plot(squeeze(std(cluster_out)))
title('GFP of clusters on PCA space')
subplot(222)

boxplot(response,id)
p = ranksum(response(id==1),response(id==2));
title(['Response, p = ' num2str(p,1)])
subplot(223)
imagesc(cluster_out(:,:,1),[-4 4])
subplot(224)
imagesc(cluster_out(:,:,2),[-4 4])


%%% reconstruct clusters to the sensor space 

for i = 1:size(cluster_out,3)
    
   
    cluster_reconstructed(:,:,i) = coeff(:,n_comp)*cluster_out(:,:,i);
    
    
end

figure('Name','All Intensities - GFP of clusters on reconstructed space','NumberTitle','off');

subplot(221)
plot(squeeze(std(cluster_reconstructed)))
title('GFP of clusters on reconstructed space')
subplot(222)
boxplot(intensity,id)

p = ranksum(intensity(id==1),intensity(id==2));
title(['Intensity, p = ' num2str(p,1)])

subplot(223)
imagesc(cluster_reconstructed(:,:,1),[-1.2 1.2])
subplot(224)
imagesc(cluster_reconstructed(:,:,2),[-1.2 1.2])


%% middle intensity 3
%%%% repeat only for middle intensity 3 

data_pca_mid = data_pca(:,:,intensity==mid_intensity );
id_mid = id(intensity==mid_intensity);
response_mid = response(intensity==mid_intensity);
intensity_mid = intensity(intensity==mid_intensity);

for i = 1:n_cluster

cluster_out_mid(:,:,i) = mean(data_pca_mid(:,:,id_mid==i),3);

end


for i = 1:n_cluster
   
    cluster_reconstructed_mid(:,:,i) = coeff(:,n_comp)*cluster_out_mid(:,:,i);
    
    
end


figure('Name','Middle Intensity - GFP of clusters on PCA space','NumberTitle','off');

subplot(221)
plot(squeeze(std(cluster_out_mid)))
title('GFP of clusters on PCA space')
subplot(222)
boxplot(response_mid,id_mid)
p = ranksum(response_mid(id_mid==1),response_mid(id_mid==2));
title(['Response, p = ' num2str(p,1)])
subplot(223)
imagesc(cluster_out_mid(:,:,1),[-4 4])
subplot(224)
imagesc(cluster_out_mid(:,:,2),[-4 4])




figure('Name','Middle Intensity - GFP of clusters on reconstructed space','NumberTitle','off');

subplot(221)
plot(squeeze(std(cluster_reconstructed_mid)))
title('GFP of clusters on reconstructed space')
subplot(222)
boxplot(intensity_mid,id_mid)
title('Intensity')
subplot(223)
imagesc(cluster_reconstructed_mid(:,:,1),[-1.2 1.2])
subplot(224)
imagesc(cluster_reconstructed_mid(:,:,2),[-1.2 1.2])

%% %%%% repeat only for middle intensities 3-4 


data_pca_mid = data_pca(:,:,intensity==mid_intensity | intensity==mid_intensity+1);
id_mid = id(intensity==mid_intensity | intensity==mid_intensity+1);
response_mid = response(intensity==mid_intensity | intensity==mid_intensity+1);
intensity_mid = intensity(intensity==mid_intensity | intensity==mid_intensity+1);

for i = 1:n_cluster

cluster_out_mid(:,:,i) = mean(data_pca_mid(:,:,id_mid==i),3);

end


for i = 1:n_cluster
   
    cluster_reconstructed_mid(:,:,i) = coeff(:,n_comp)*cluster_out_mid(:,:,i);
    
    
end


figure('Name','Middle Intensities - GFP of clusters on PCA space','NumberTitle','off');

subplot(221)
plot(squeeze(std(cluster_out_mid)))
title('GFP of clusters on PCA space')
subplot(222)
boxplot(response_mid,id_mid)
p = ranksum(response_mid(id_mid==1),response_mid(id_mid==2));
title(['Response, p = ' num2str(p,1)])
subplot(223)
imagesc(cluster_out_mid(:,:,1),[-4 4])
subplot(224)
imagesc(cluster_out_mid(:,:,2),[-4 4])




figure('Name','Middle Intensities - GFP of clusters on reconstructed space','NumberTitle','off');

subplot(221)
plot(squeeze(std(cluster_reconstructed_mid)))
title('GFP of clusters on reconstructed space')
subplot(222)
boxplot(intensity_mid,id_mid)
p = ranksum(intensity_mid(id_mid==1),intensity_mid(id_mid==2));
title(['Intensity, p = ' num2str(p,1) ', means ' num2str(mean(intensity_mid(id_mid==1))) '  ' num2str(mean(intensity_mid(id_mid==2)))])
subplot(223)
imagesc(cluster_reconstructed_mid(:,:,1),[-1.2 1.2])
subplot(224)
imagesc(cluster_reconstructed_mid(:,:,2),[-1.2 1.2])





