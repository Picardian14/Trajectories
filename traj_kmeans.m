function [id_out, clust_out, distortion, silhouette_all, c_out, clust_size_out] = traj_kmeans(data,n_clust,nrep)
% data structure = dimensions x time x channels

min_dist = inf;

seed_method = '++';
%seed_method = 'rand_mean';


parfor n = 1:nrep

%%% seed 
id(n,:) = randi(n_clust,size(data,3),1);
clust = nan(size(data,1),size(data,2),n_clust);
dist = nan(size(data,3),n_clust);
dist_init = nan(size(data,3),n_clust);

%%% Different methods for selecting initial cluster
switch seed_method
    case 'rand_mean'     
        for i = 1:n_clust
            clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),20,'round',3); %% mean/median/trm
        end
    case 'rand'
        
    case '++'
        %%% initial cluster
        %disp('seeding')     
        index = randi(size(data,3),1);
        % disp(['index 1: ' num2str(index)])
        clust(:,:,1) = data(:,:,index);       
        for clust_index = 1:(n_clust-1)
            mat = data - repmat(clust(:,:,clust_index),[1 1 size(data,3)]);
            dist_init(:,clust_index) = squeeze(sum(sum(abs(mat))));
            proba = cumsum((min(dist_init')).^2)./sum((min(dist_init')).^2);
            index = find(proba>rand);
            index = index(1);
            %disp(['index ' num2str(clust_index+1) ': ' num2str(index)])
            clust(:,:,clust_index+1) = data(:,:,index);
        end
end


aux = zeros(1,size(data,3));

%%% trajectory k-means 
for rep = 1:2000    
    for i=1:n_clust
        mat = data - repmat(clust(:,:,i),[1 1 size(data,3)]);
        dist(:,i) = squeeze(sum(sum(abs(mat))));
    end 
    [a,id(n,:)] = min(dist'); 
    %id(n, :) = randi(n_clust, size(id(n, :)));
   % disp(['iter: ' num2str(rep) ' sum dist: ' num2str(sum(a))])
    
    if prod(id(n,:)==aux)==1
        break
    end
    for i = 1:n_clust
        clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),20,'round',3); %% mean/median/trm
        %sum_d(i) = sum(data(:, :, id(n, :)==i) - repmat(clust(:,:,i),[1 1 size(data,3)]))
    end
    
    aux = id(n,:);
       
end



sum_dist(n) = sum(a);

tmpvarc = nan(size(data,3),size(data,3),1);
% calculate connectivity matrix
for i_elem = 1:size(aux,2)
    for j_elem = 1:size(aux,2)
        if aux(i_elem) == aux(j_elem)
            tmpvarc(i_elem,j_elem) = 1;
        else
            tmpvarc(i_elem,j_elem) = 0;
        end
    end
end

c(:, :, n) = tmpvarc;

%disp(['rep: ' num2str(n) ', iter: ' num2str(rep) ', total distance:' num2str(sum_dist(n))])
    
end


[a, n] = min(sum_dist);
clust_sizes=zeros(n_clust, nrep);
disp(['Distance min:' num2str(a)])
for i_clust=1:n_clust
    for i_rep=1:nrep
        clust_sizes(i_clust,i_rep) = sum(id(i_rep, :)==i_clust);
    end
end
clust_size_out = clust_sizes;

id_out = squeeze(id(n,:));

for i = 1:n_clust
    clust_out(:,:,i) = trimmean(data(:,:,id_out==i),20, 'round', 3); %% mean/median/trm    
end
distortion = 0;
for i = 1:n_clust
    mat = data(:,:,id_out==i) - repmat(clust_out(:,:,i),[1 1 size(data(:,:,id_out==i),3)]);
    distortion = distortion + squeeze(sum(sum(sum(abs(mat)))));
end

disp("Silhouette for "+ n_clust+ " clusters");
if n_clust==2
    n_obs = size(data(:, :, id_out==1),3);
    observations = data(:, :, id_out==1);
    tmp_clust = clust_out(:, :, 2);
    for j_obs = 1:n_obs
            sample = observations(:, :, j_obs);
            % intra cluster distance
            % Mean over the number of observations in the cluster
            manhattandist = mean(observations - repmat(sample, [1 1 size(observations, 3)]), 3);
            a = sum(sum(abs(manhattandist)));
            % nearest cluster distance. Did With centroid               
            b = sum(sum(abs(sample - tmp_clust)));        
            % a and b are distances, so i just use Manhattan distance for all
            % time points
            silhouette_clust(:, j_obs) = (b-a)/(max(a,b));            
    end
    disp("Clust " + 1 + " score: " + mean(silhouette_clust)+ "// #Elements: "+n_obs);
    silhouette_all(:, 1) = mean(silhouette_clust);

    n_obs = size(data(:, :, id_out==2),3);
    observations = data(:, :, id_out==2);
    tmp_clust = clust_out(:, :, 1);
    for j_obs = 1:n_obs
            sample = observations(:, :, j_obs);
            % intra cluster distance
            % Mean over the number of observations in the cluster
            manhattandist = mean(observations - repmat(sample, [1 1 size(observations, 3)]), 3);
            a = sum(sum(abs(manhattandist)));
            % nearest cluster distance. Did With centroid               
            b = sum(sum(abs(sample - tmp_clust)));        
            % a and b are distances, so i just use Manhattan distance for all
            % time points
            silhouette_clust(:, j_obs) = (b-a)/(max(a,b));            
    end
    disp("Clust " + 2 + " score: " + mean(silhouette_clust)+ "// #Elements: "+n_obs);
    silhouette_all(:, 2) = mean(silhouette_clust);    
else
    for i_clust = 1:n_clust
        n_obs = size(data(:, :, id_out==i_clust),3);
        observations = data(:, :, id_out==i_clust);
        % Get nearest clusters for current clusters
        I = [ i_clust ];
        tmp_clust = clust_out(:, :, setdiff(1:end, I));
        % 1. Get nearest for all samples in the current cluster of observations
        for tmp_i=1:n_clust-1        
            mat = observations - repmat(tmp_clust(:,:,tmp_i),[1 1 size(observations,3)]);
            dist(:,tmp_i) = squeeze(sum(sum(abs(mat))));
        end    
        % The index of the nearest centroid for each observation in the current
        % cluster is in tmp_id
        [a,tmp_id] = min(dist');     
        
        for j_obs = 1:n_obs
            sample = observations(:, :, j_obs);
            % intra cluster distance
            % Mean over the number of observations in the cluster
            manhattandist = mean(observations - repmat(sample, [1 1 size(observations, 3)]), 3);
            a = sum(sum(abs(manhattandist)));
            % nearest cluster distance. Did With centroid               
            b = sum(sum(abs(sample - tmp_clust(:, :, tmp_id(j_obs)))));        
            % a and b are distances, so i just use Manhattan distance for all
            % time points
            silhouette_clust(:, j_obs) = (b-a)/(max(a,b));
        end
        disp("Clust " + i_clust + " score: " + mean(silhouette_clust)+ "// #Elements: "+n_obs);
        silhouette_all(:, i_clust) = mean(silhouette_clust);
        clear dist;
    end
end
disp("Clust Overall score: " + mean(silhouette_all));
c_out = c;
end


