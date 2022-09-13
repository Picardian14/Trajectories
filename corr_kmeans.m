function [id_out, clust_out] = topo_kmeans(data,n_clust,nrep)
% data structure = dimensions x time x channels

min_dist = inf;

seed_method = '++';
%seed_method = 'rand_mean';

verb = 1; 

parfor n = 1:nrep

%%% seed 
id(n,:) = randi(n_clust,size(data,3),1);

clust = nan(size(data,1),size(data,2),n_clust);
dist = nan(size(data,3),n_clust);
dist_init = nan(size(data,3),n_clust);
mat = nan(size(data,1),size(data,3));
switch seed_method
    case 'rand_mean'

for i = 1:n_clust
    clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),20,'round',3); %% mean/median/trm
end

    case 'rand'
        
    case '++'
        %%% initial cluster
if verb; disp('seeding'); end
     
index = randi(size(data,3),1);
if verb;     disp(['index 1: ' num2str(index)]); end
     clust(:,:,1) = data(:,:,index);
   
     for clust_index = 1:(n_clust-1)
      aux = clust(:,:,clust_index);
%         mat = data - repmat(aux,[1 1 size(data,3)]);

X = permute(data,[2 1 3]);
Z = corr(aux', X(:,:));

for jjj = 1:size(data,1)
aux_d =  (ones(size(data,3),1) - Z(jjj,jjj:4:end)')/2;
%mat(jjj,:) =  0.5*log((1+aux_d)./(1-aux_d));
mat(jjj,:) =  aux_d;

end


      dist_init(:,clust_index) = squeeze(sum(mat));
      proba = cumsum((min(dist_init')).^2)./sum((min(dist_init')).^2);
      index = find(proba>rand);
      index = index(1);
if verb;           disp(['index ' num2str(clust_index+1) ': ' num2str(index)]); end
      clust(:,:,clust_index+1) = data(:,:,index);
     end
end

if verb; disp('kmeans'); end

aux_id = zeros(1,size(data,3));

%%% trajectory k-means 
for rep = 1:500
    
    for i=1:n_clust
   % mat = data - repmat(clust(:,:,i),[1 1 size(data,3)]);
      aux = clust(:,:,i);

X = permute(data,[2 1 3]);
Z = corr(aux', X(:,:));

for jjj = 1:size(data,1)
aux_d =  (ones(size(data,3),1) - Z(jjj,jjj:4:end)')/2;
%mat(jjj,:) =  0.5*log((1+aux_d)./(1-aux_d));
mat(jjj,:) =  aux_d;

end

    dist(:,i) = squeeze(sum(mat));
    end
 
    [a,id(n,:)] = min(dist');
    
    
    
if 1;    disp(['rep: ' num2str(rep) ' sum dist: ' num2str(sum(a))]); end

    if sum(abs(id(n,:)-aux_id))==0
        break
    end
    
for i = 1:n_clust
    clust(:,:,i) = trimmean(data(:,:,id(n,:)==i),20,'round',3); %% mean/median/trm
end

    aux_id = id(n,:);
   
end

sum_dist(n) = sum(a);

disp(['rep: ' num2str(n) ', iter: ' num2str(rep) ', total distance:' num2str(sum_dist(n))])
    
end


[a, n] = min(sum_dist);

disp(['Distance min:' num2str(a)])
 

id_out = squeeze(id(n,:));

for i = 1:n_clust
    clust_out(:,:,i) = trimmean(data(:,:,id_out==i),20, 'round', 3); %% mean/median/trm
end



end


