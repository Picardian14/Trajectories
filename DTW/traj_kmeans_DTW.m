function [id_out, clust_out,a] = traj_kmeans(data,n_clust,nrep,limit,thr)
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

switch seed_method
    case 'rand_mean'

for i = 1:n_clust
    clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),5,'round',3); %% mean/median/trm
end

    case 'rand'
        
    case '++'
        %%% initial cluster
%disp('seeding')
     
index = randi(size(data,3),1);
    % disp(['index 1: ' num2str(index)])
     clust(:,:,1) = data(:,:,index);
   
     for clust_index = 1:(n_clust-1)
%         mat = movmean(data,m_avg,2,'omitnan') - repmat(clust(:,:,clust_index),[1 1 size(data,3)]);
        

mat = data - repmat(clust(:,:,clust_index),[1 1 size(data,3)]);
      
      dist_init(:,clust_index) = squeeze(sum(sum(abs(mat))));         
         

proba = cumsum((min(dist_init')).^2)./nansum((min(dist_init')).^2);
      index = find(proba>rand);
      index = index(1);
     %      disp(['index ' num2str(clust_index+1) ': ' num2str(index)])
      clust(:,:,clust_index+1) = data(:,:,index);
     end
end


aux = zeros(1,size(data,3));

%%% trajectory k-means 
for rep = 1:250
    
    for i=1:n_clust
  
    
  for trial = 1:size(data,3)
      
      dist(trial,i) = dtw(data(:,:,trial),clust(:,:,i),limit,'absolute');
      
  end

    end

 
    [a,id(n,:)] = min(dist');
    
    
    disp(['iter: ' num2str(rep) ' sum dist: ' num2str(sum(a))])
    
    if prod(id(n,:)==aux)==1
        break
    end
    
for i = 1:n_clust
    clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),thr,'round',3); %% mean/median/trm
end



    aux = id(n,:);
   
end

sum_dist(n) = sum(a);

disp(['-- rep: ' num2str(n) ', iter: ' num2str(rep) ', total distance:' num2str(sum_dist(n))])
    
end


[a, n] = min(sum_dist);

disp(['------- Distance min:' num2str(a)])
 

id_out = squeeze(id(n,:));

for i = 1:n_clust
    clust_out(:,:,i) = trimmean(data(:,:,id_out==i),thr, 'round', 3); %% mean/median/trm
end



end


