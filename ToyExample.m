%% Setup points
clear all
close all


X11 = 1;
Y11 = 1;

X12 = 6;
Y12 = 1;

X21 = 3;
Y21 = 0;

X22 = 3;
Y22 = 7;
% Two dimensions or sensors?

X31 = 1.5;
Y31 = 1;

N = 100; % Number or trials or trajectories

A(:, 1) = X11 + rand(N,1)/4;
A(:, 2) = Y11 + rand(N,1)/4;

B(:, 1) = X12 + rand(N,1)/4;
B(:, 2) = Y12 + rand(N,1)/4;

Z(:, 1) = X21 + rand(N,1)/4;
Z(:, 2) = Y21 + rand(N,1)/4;

W(:, 1) = X22 + rand(N,1)/4;
W(:, 2) = Y22 + rand(N,1)/4;

V(:, 1) = X31 + rand(N,1)/4;
V(:, 2) = Y31 + rand(N,1)/4;

%% PLOT points
figure
hold on
plot(A(1,1), A(1,2), 'or', Color='blue');
plot(B(1,1), B(1,2), 'or', Color='red');
plot(Z(1,1), Z(1,2), 'or', Color='green');
plot(W(1,1), W(1,2), 'or', Color='black');
plot(V(1,1), V(1,2), 'or', Color='cyan');
hold off

Tmax = 150; % Timepoints
%% Generate
for i = 1:N

    ABs(:, i, 1) = linspace(A(i,1) + rand(1,1)/4 * randi([-1 1]), B(i, 1)+ rand(1,1)/4 * randi([-1 1]), Tmax); % X
    ABs(:, i, 2) = linspace(A(i,2) + rand(1,1)/4 * randi([-1 1]), B(i, 2)+ rand(1,1)/4 * randi([-1 1]), Tmax); % Y
    ZWs(:, i, 1) = linspace(Z(i,1) + rand(1,1)/4 * randi([-1 1]), W(i, 1)+ rand(1,1)/4 * randi([-1 1]), Tmax); % X
    ZWs(:, i, 2) = linspace(Z(i,2) + rand(1,1)/4 * randi([-1 1]), W(i, 2)+ rand(1,1)/4 * randi([-1 1]), Tmax); % Y
end

%% Plot
figure
hold on
for j = 1:N
    plot(ABs(:, j, 1), ABs(:, j, 2), 'color', 'blue');
    plot(ZWs(:, j, 1), ZWs(:, j, 2), 'color', 'black');
end
hold off

%% Cluster
data = cat(2, ZWs, ABs);
data = permute(data, [3,1,2]); % Dimensions % time % channels
n_clust = 2;
n_rep = 2;
[id, clust] = traj_kmeans(data, n_clust, n_rep);

%% Plot cluters
figure
hold on
plot(clust(1, :, 1),clust(2, :, 1),'color', 'red')
plot(clust(1, :, 2),clust(2, :, 2),'color', 'black')

%%xlim([0,8])
%%ylim([0 8])
hold off

%% More complex data

% WA changes at half of the way
idxaux = 1:N;
idxlist1 = idxaux(randperm(N));
idxlist2 = idxaux(randperm(N));
for i = 1:N
    Ztemp(i, 1) = Z(idxlist1(i), 1) + rand(1,1)/4 * randi([-1 1]);
    Ztemp(i, 2) = Z(idxlist1(i), 2) + rand(1,1)/4 * randi([-1 1]);
    Wtemp(i, 1) = W(idxlist1(i), 1) + rand(1,1)/4 * randi([-1 1]);
    Wtemp(i, 2) = W(idxlist1(i), 2) + rand(1,1)/4 * randi([-1 1]);
    tempWAs(:, i, 1) = linspace(Wtemp(i,1), Ztemp(i, 1), Tmax/2); % X
    tempWAs(:, i, 2) = linspace(Wtemp(i,2), Wtemp(i, 2)/2, Tmax/2); % Y    
end

for i = 1:N
    WAs(:, i, 1) = cat(1,tempWAs(:, i, 1),linspace(Ztemp(i,1), A(idxlist1(i), 1), Tmax/2)'); % X
    WAs(:, i, 2) = cat(1,tempWAs(:, i, 2),linspace(Wtemp(i,2)/2, A(idxlist1(i), 2), Tmax/2)'); % Y    
end

%V is another middle point, closer to A
%I will make a WV point 
clear Ztemp;
for i = 1:N
    Ztemp(i, 1) = Z(idxlist2(i), 1) + rand(1,1)/4 * randi([-1 1]);
    Ztemp(i, 2) = Z(idxlist2(i), 2) + rand(1,1)/4 * randi([-1 1]);
    tempWVs(:, i, 1) = linspace(W(idxlist2(i),1)+ rand(1,1)/4 * randi([-1 1]), Ztemp(i, 1), Tmax/2); % X
    tempWVs(:, i, 2) = linspace(W(idxlist2(i),2), W(idxlist2(i), 2)/2, Tmax/2); % Y    
end

for i = 1:N
    WVs(:, i, 1) = cat(1,tempWVs(:, i, 1),linspace(Ztemp(i,1), V(idxlist2(i), 1), Tmax/2)'); % X
    WVs(:, i, 2) = cat(1,tempWVs(:, i, 2),linspace(W(idxlist2(i),2)/2, V(idxlist2(i), 2), Tmax/2)'); % Y    
end

%% Plot complex data

figure
hold on
for j = 1:N
    plot(WVs(:, j, 1), WVs(:, j, 2), 'color', 'blue');
    plot(WAs(:, j, 1), WAs(:, j, 2), 'color', 'black');
end
hold off
%% Cluster
data = cat(2, WAs, WVs);
data = permute(data, [3,1,2]); % Dimensions % time % channels
n_clust = 2;
[id, clust] = traj_kmeans(data, 2, 2);

%% Plot clusters
figure
hold on
plot(clust(1, :, 1),clust(2, :, 1),'color', 'red')
plot(clust(1, :, 2),clust(2, :, 2),'color', 'black') 

%%xlim([0,8])
%%ylim([0 8])
hold off


%% Evaluation
% Cluster 2, black should go with WV lines
% Cluster 1, red should go with WA lines
for i=1:n_clust
    mat = WVs - permute(repmat(clust(:,:,i),[1 1 size(WVs,2)]), [2,3,1]);
    dist(:,i) = squeeze(sum(sum(abs(permute(mat, [3 1 2])))));
end
%%
[a_wv,id_wv(1,:)] = min(dist');
%sum(id_wv - 1)/100
%97.5
%%
for i=1:n_clust
    mat_wa = WAs - permute(repmat(clust(:,:,i),[1 1 size(WAs,2)]), [2,3,1]);
    dist_wa(:,i) = squeeze(sum(sum(abs(permute(mat_wa, [3 1 2])))));
end
%%
[a_wa,id_wa(1,:)] = min(dist_wa');
%1 - sum(id_wa - 1)/100
%95
%% Shuffled

%data_shuf = data(randperm(size(data,2)), :, :);
%[id_shuf, clust_shuf] = traj_kmeans(data_shuf, 2, 2);

%% Time divergent

clear WAs
clear WVs
clear tempWVs
clear tempWAs
% WA changes at half of the way
idxaux = 1:N;
idxlist1 = idxaux(randperm(N));
idxlist2 = idxaux(randperm(N));
for i = 1:N
    WZs(:, i, 1) = linspace(W(idxlist1(i),1), Z(idxlist1(i), 1), Tmax); % X
    WZs(:, i, 2) = linspace(W(idxlist1(i),2), Z(idxlist1(i), 2), Tmax); % Y    
end

% Second trajectory will just be more spaced, and repeat the end point
% We'll see if it makes sense

for i = 1:N
    WX = W(idxlist2(i),1);
    WY = W(idxlist2(i),2);
    tempWAfast(:, i, 1) = linspace(WX, (Z(idxlist2(i), 1)-WX)/2, Tmax/2); % X
    tempWAfast(:, i, 2) = linspace(WY, WY/2, Tmax/2); % Y    
    middle_point(i, 1) = (Z(idxlist2(i), 1)-WX)/2;
    middle_point(i, 2) = WY/2;
    tempWAfast2(:, i, 1) = linspace(middle_point(i, 1), Z(idxlist2(i), 2), round(Tmax/4));
    tempWAfast2(:, i, 2) = linspace(middle_point(i, 2), Z(idxlist2(i), 2), round(Tmax/4));
    end_point(:, i, 1) = repmat(Z(idxlist2(i), 1), 1,Tmax-(round(Tmax/4)+Tmax/2));
    end_point(:, i, 2) = repmat(Z(idxlist2(i), 2), 1,Tmax-(round(Tmax/4)+Tmax/2));
end

%%
for i = 1:N    
    WZfast(:, i, 1) = cat(1,tempWAfast(:, i, 1),tempWAfast2(:, i, 1), end_point(:, i, 1)); % X
    WZfast(:, i, 2) = cat(1,tempWAfast(:, i, 2),tempWAfast2(:, i, 2), end_point(:, i, 2)); % X
end

%% Plot time data
% DOING IT WRONG
figure
hold on
for j = 1:N
    plot(WZs(:, j, 1), WZs(:, j, 2), 'color', 'blue');
    plot(WZfast(:, j, 1), WZfast(:, j, 2), 'color', 'black');
end
hold off




