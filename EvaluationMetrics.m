%% Elbox method
% 1. Get the sum of distances between all points in a cluster and the centroid
% 2. Save the number of clusters and the sum of distances, for each k
% 3. Plot k x distance

% TODO: DO IT RIGHT. ES LA SUMA DE TODAS LAS DISTANCIAS AL FINAL? NO TIENE
% SENTUDO PORQUE TE QUEDA EL QUE TIENE MAS CLUSTERS COMO EL DE MAS
% DISTANCIA


n_rep = 4;
n_clust = 16;
start_clust = 3;
for k=start_clust:n_clust
    [id, clust, distortion] = traj_kmeans(meandata, k, n_rep);
    
    for i = 1:k
        mat_out = meandata - repmat(clust(:,:,i),[1 1 size(meandata,3)]);
        dist_out(:, i) = squeeze(sum(sum(sum(abs(mat_out)))));
    end
    D(1, k) = k;
    D(2, k) = distortion;%sum(dist_out);
end
%% Plot
% ELBOW IS NOT VERY ELBOW. MAYBE TIME FACTOR IS NOT ELBOW COMPATIBLE
figure
hold on;
plot(D(1,start_clust:end),D(2,start_clust:end));    
plot(D(1,start_clust:end),D(2,start_clust:end), 'or');
hold off

%%

figure
tiledlayout(1,2)
nexttile
colors = ["red", "black", "green", "blue"];
hold on
for i_clust = 1:n_clust
    plot(clust(1, :, i_clust),'color', colors(1,i_clust));                        
end
title('Subject '+string(roi+pad))
hold off     
nexttile

hold on
for i_clust = 1:n_clust
    plot(clust(2, :, i_clust),'color', colors(1,i_clust));                        
end
title('Subject '+string(roi+pad))
hold off

parcelsxfig = size(clust,1)/2;
for ifig = 1:2
    pad = (ifig-1)*parcelsxfig;
    figure
    tiledlayout(5,3);    
    colors = ["red", "black", "green", "blue"];
    for roi = 1:parcelsxfig
        nexttile;
        hold on
        for i_clust = 1:n_clust
            plot(clust(roi+pad, :, i_clust),'color', colors(1,i_clust));                        
        end
        title('Subject '+string(roi+pad))
        hold off     
    end
end

%% EXPLORING K
data_reshape = tsdata_pca;
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
ks = start_clust:n_clust
plot(ks,rho(start_clust:end));    
plot(ks,rho(start_clust:end), 'or');
hold off
title("Dispersion")