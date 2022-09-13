clear all
close all


load allp_mean_env_Valid_SSOAipsi_contra_LSOAipsi_contra_Invalid_SSOAipsi_contra_LSOAipsi_contra_cleaneddata5cons.mat

n_clusters = 6;


me = mean(allp_mean_env(:,:),2);
st = std(allp_mean_env(:,:),[],2);


for j = 1:size(allp_mean_env,1)
    
    for i = 1:size(allp_mean_env,3)
        
   allp_mean_env(j,:,i) = (allp_mean_env(j,:,i) - me(j))/st(j);     
        
    end
    
end


data = permute(allp_mean_env, [3 2 1]);

dims = size(data,1);
z = 0;

for DTW = [1 2 3 4 5 6 8 10 12 15 20 30 40 50 75 100]

    z = z +1;
[trials_id,clust, dist] = traj_kmeans_DTW(data,n_clusters,500,DTW,10);

save(['DTW' num2str(DTW) '.png'],'trials_id','clust','dist')


clust_all(:,:,:,z) = clust;

for i = 1:n_clusters
    
    n(i) = sum(trials_id==i);
end

[a,b] = sort(n);


n = 0;
for d = 1:dims
for i = 1:n_clusters
    n =n +1;
    subplot(dims,n_clusters,n)
    
    plot(squeeze(clust(d,:,b(i))))
    
    if d ==1
        title(['# = ' num2str(sum(trials_id==b(i)))])
    end
    axis([0 114 -1 2.5])
end
end

saveas(gcf,['DTW' num2str(DTW) '.png'])
close

end

clust_order = clust_all;

DTW_1 = 1;
    for DTW_2 = (DTW_1+1):z
        
        d = nan(n_clusters);
for i_1 = 1:n_clusters
for i_2 = 1:n_clusters

X = clust_all(:,:,i_1,DTW_1);
Y = clust_all(:,:,i_2,DTW_2);
d(i_1,i_2) = sum(abs(X(:)-Y(:)));     
                
    end
end

for i = 1:n_clusters
    
[a,b] = min(d(:));
[x,y] = ind2sub([n_clusters,n_clusters],b);
map(x,DTW_2-1) = y ;
d(x,:) = nan;
d(:,y) = nan;

clust_order(:,:,x,DTW_2) = clust_order(:,:,y,DTW_2);


end

end

    
   
    
    
n = 0;
for d = 1:dims
for i = 1:n_clusters
    n =n +1;
    subplot(dims,n_clusters,n)
    
    
    for DTW = 1:3
    
    plot(squeeze(clust_order(d,:,i,DTW)))
    
    hold on
    
    end
    if d ==1
    legend
    end
    axis([0 114 -1 2.5])
end
end




