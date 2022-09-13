%% Setup points
clear all
close all

hard = false;
if hard
    X11 = 1.25;
    Y11 = 0.75;
    
    X12 = 1.3;
    Y12 = 0.7;
    
    X21 = 0.8;
    Y21 = 0.7;
    
    X22 = 1;
    Y22 = 0.2;
    
    X31 = 1.5;
    Y31 = 1;
else
    % Easier config
    X11 = 1;
    Y11 = 0.85;
    
    X12 = 1.3;
    Y12 = 0.5;
    
    X21 = 0.8;
    Y21 = 0.6;
    
    X22 = 1.1;
    Y22 = 0.1;
    
    X31 = 1.5;
    Y31 = 1;
end

N = 100; % Number or trials or trajectories

A(:, 1) = X11 + rand(N,1)/3;
A(:, 2) = Y11 + rand(N,1)/3;

B(:, 1) = X12 + rand(N,1)/3;
B(:, 2) = Y12 + rand(N,1)/3;

Z(:, 1) = X21 + rand(N,1)/3;
Z(:, 2) = Y21 + rand(N,1)/3;

W(:, 1) = X22 + rand(N,1)/3;
W(:, 2) = Y22 + rand(N,1)/3;

V(:, 1) = X31 + rand(N,1)/3;
V(:, 2) = Y31 + rand(N,1)/3;

data_t1(:, 1) = cat(1,cat(1,cat(1, A(1:50, 1),B(1:50, 1)),W(1:50, 1)), V(1:50, 1));
data_t1(:, 2) = cat(1,cat(1,cat(1, A(1:50, 2),B(1:50, 2)),W(1:50, 2)), V(1:50, 2));
data_t2(:, 1) = cat(1,cat(1,cat(1, A(51:100, 1),B(51:100, 1)),W(51:100, 1)), V(51:100, 1));
data_t2(:, 2) = cat(1,cat(1,cat(1, A(51:100, 2),B(51:100, 2)),W(51:100, 2)), V(51:100, 2));
%data = cat(1, data_t1, data_t2);
data(:, 1, :) = data_t1;
data(:, 2, :) = data_t2;
zdata = zscore(data, 0, 1);
data = permute(data, [3 2 1]);
zdata = permute(zdata, [3 2 1]);
% zdata is timexdimsxobs
figure
hold on
scatter(A(:,1), A(:,2), [], 'red');
scatter(B(:,1), B(:,2), [], 'blue');
scatter(W(:,1), W(:,2), [], 'green');
scatter(V(:,1), V(:,2), [], 'black');
title('original')

hold off
%%
figure
hold on
% zdata is timexdimsxobs
scatter(squeeze(zdata(1,1,1:50)), squeeze(zdata(1,2,1:50)), [], 'red');
scatter(squeeze(zdata(1,1,51:100)), squeeze(zdata(1,2,51:100)), [], 'blue');
scatter(squeeze(zdata(1,1,101:150)), squeeze(zdata(1,2,101:150)), [], 'green');
scatter(squeeze(zdata(1,1,151:200)), squeeze(zdata(1,2,151:200)), [], 'black');
title('zscored time 1')
hold off
%%
figure
hold on
scatter(squeeze(zdata(2,1,1:50)), squeeze(zdata(2,2,1:50)), [], 'red');
scatter(squeeze(zdata(2,1,51:100)), squeeze(zdata(2,2,51:100)), [], 'blue');
scatter(squeeze(zdata(2,1,101:150)), squeeze(zdata(2,2,101:150)), [], 'green');
scatter(squeeze(zdata(2,1,151:200)), squeeze(zdata(2,2,151:200)), [], 'black');
title('zscored time 2')
hold off



%% Clustering 4 groups with mixed
n_rep = 128;
n_clust = 8;
clusts = {};
best_sil = 0;
best_clust = {};
clust_obs={};
clust_size={};
% colors = ["red", "blue", "green", "black"];
for i_clust=2:n_clust
    cur_clust = {};
    [id, clust, distortion, silhouette, c, clust_sizes] = traj_kmeans(zdata, i_clust, n_rep);
    if best_sil<mean(silhouette)            
        best_clust = {id, clust, silhouette};
        best_sil=mean(silhouette);               
    end
    C = mean(c,3);
    % calculate dispersion coefficient as in (Kim & Park, Bioinformatics, 2007)
    tmp = 0;
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            tmp = tmp + (4 * (C(i,j) - 1/2)^2);
        end
    end
    rho(i_clust) = 1/((size(C,2))^2) * tmp;
    figure
    hold on
    leg_text = "N_clust: ";
    colors = [ "black", "white", "red", "green", "blue", "cyan", "magenta", "yellow", "gray", "lightBlue", "orange", "darkGreen"];
    for j_clust=1:i_clust
        cl = zdata(:, :, id==j_clust);        
        cur_clust{end+1} = {cl, size(cl, 3)};        
        scatter(squeeze(cl(1,1,:)), squeeze(cl(1,2,:)), [], colors(j_clust));        
        leg_text=leg_text + cur_clust{end}{2} + ", ";
    end
    leg_text=leg_text + "Rho: " + rho(i_clust)
    title(leg_text)
    hold off
    clust_obs{end+1} = cur_clust;  
    clust_size{end+1} = clust_sizes;
end

%%

for t=1:size(clust_size, 2)
    cur_clusts = clust_size{t};
    figure
    tiledlayout(1, size(cur_clusts, 1));
    for l=1:size(cur_clusts, 1)
        nexttile;
        histogram(cur_clusts(l, :));
    end
end

for t=1:size(clust_size, 2)
    sum(std(clust_size{t}'))      
end
%figure
%hold on 
%for i_clust=1:n_clust
%    plot(linspace(1, 16, 16), squeeze(nclusts(:, i_clust)), 'color', colors(i_clust))
%end
%hold off

%% Preparing two easy
clear data_t1;
clear data_t2;
data_t1(:, 1) = cat(1,cat(1,cat(1, W(1:50,1),V(1:50,1))));
data_t1(:, 2) = cat(1,cat(1,cat(1, W(1:50,2),V(1:50,2))));
data_t2(:, 1) = cat(1,cat(1,cat(1, W(51:100,1),V(51:100,1))));
data_t2(:, 2) = cat(1,cat(1,cat(1, W(51:100,2),V(51:100,2))));
%data = cat(1, data_t1, data_t2);
clear WVdata
clear WVzdata
WVdata(:, 1, :) = data_t1;
WVdata(:, 2, :) = data_t2;
WVzdata = zscore(WVdata, 0, 1);
WVdata = permute(WVdata, [3 2 1]);
WVzdata = permute(WVzdata, [3 2 1]);

figure
hold on
% zdata is timexdimsxobs
scatter(squeeze(WVzdata(1,1,1:50)), squeeze(WVzdata(1,2,1:50)), [], 'red');
scatter(squeeze(WVzdata(1,1,51:100)), squeeze(WVzdata(1,2,51:100)), [], 'blue');
title('zscored time 1')
hold off
figure
hold on
scatter(squeeze(WVzdata(2,1,1:50)), squeeze(WVzdata(2,2,1:50)), [], 'red');
scatter(squeeze(WVzdata(2,1,51:100)), squeeze(WVzdata(2,2,51:100)), [], 'blue');
title('zscored time 2')
hold off

%% Clustering two easy groups
[id, clust, distortion, eval, c, clust_sizes] = traj_kmeans(WVzdata, 2, 64);
c1 = WVzdata(:, :, id==1);
c2 = WVzdata(:, :, id==2);       
figure
hold on
% zdata is timexdimsxobs
scatter(squeeze(c1(1,1,:)), squeeze(c1(1,2,:)), [], 'red');
scatter(squeeze(c2(1,1,:)), squeeze(c2(1,2,:)), [], 'blue');    
title("Time 1");
hold off

figure
hold on
% zdata is timexdimsxobs
scatter(squeeze(c1(2,1,:)), squeeze(c1(2,2,:)), [], 'red');
scatter(squeeze(c2(2,1,:)), squeeze(c2(2,2,:)), [], 'blue');    
title("Time 2");
hold off
%% Preparing two mixed groups
data_t1 = cat(1,cat(1,cat(1, A(1:50),B(1:50))));
data_t2 = cat(1,cat(1,cat(1, A(51:100),B(51:100))));
%data = cat(1, data_t1, data_t2);
ABdata(:, 1, :) = data_t1;
ABdata(:, 2, :) = data_t2;
ABzdata = zscore(ABdata, 0, 1);
ABdata = permute(ABdata, [3 2 1]);


figure
hold on
scatter(squeeze(ABdata(:,1,1)), squeeze(ABdata(:,2,1)), [], 'red');
scatter(squeeze(ABdata(:,1,2)), squeeze(ABdata(:,2,2)), [], 'blue');
hold off

%% Clustering two mixed groupts
[id, clust, distortion, eval] = traj_kmeans(data, 2, 1);