
for nsub=indexsub
    corr_tsdata_pca(:, :, nsub) = corrcoef(squeeze(tsdata_pca(:, :, nsub))');
end


for i=1:size(tsdata_pca, 3)
    tsdata_reconstructed(:, :, i) = coeff(:, NCOMP)*tsdata_pca(:, :, i);
end

for nsub=indexsub
    corr_tsdata_reconstructed(:, :, nsub) = corrcoef(squeeze(tsdata_reconstructed(:, :, nsub))');
end

for nsub=indexsub
    corr_tsdata(:, :, nsub) = corrcoef(squeeze(tsdata(:, :, nsub))');
end

for i=1:90
    for j=1:30
        cur_t = squeeze(tsdata(i, :, j));
        if isequal(tmp_roi_sub, cur_t)
            disp("ROI: "+i+" -- Subject: "+j);
        end
    end
end







%%


%% 


% correlate between no epochs
sample=1:2:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:4
            tmp = corrcoef(clust(i_pc, :, i_clust), clust(i_pc, :, j_clust));
            intercorr(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end


% Correlate between no epoch and epoch 100
sample=1:2:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:4
            tmp = corrcoef(clust(i_pc, sample, i_clust), clust_e(i_pc, :, j_clust));
            interepochcorr(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

%sample=1:2:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:4
            tmp = corrcoef(clust(i_pc, 101:200, i_clust), clust_e(i_pc, :, j_clust));
            interepochcorr_second(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

%correlate between no epoch and epoch 50

sample=1:4:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:2
            tmp = corrcoef(clust(i_pc, sample, i_clust), clust_e50(i_pc, :, j_clust));
            interepoch50corr(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

sample=1:4:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:2
            tmp = corrcoef(clust(i_pc, 1:50, i_clust), clust_e50(i_pc, :, j_clust));
            interepoch50corr_1(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

sample=1:4:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:2
            tmp = corrcoef(clust(i_pc, 51:100, i_clust), clust_e50(i_pc, :, j_clust));
            interepoch50corr_2(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

sample=1:4:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:2
            tmp = corrcoef(clust(i_pc, 101:150, i_clust), clust_e50(i_pc, :, j_clust));
            interepoch50corr_3(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end

sample=1:4:200;
for i_clust=1:4
    for i_pc=1:15       
        for j_clust=1:2
            tmp = corrcoef(clust(i_pc, 151:200, i_clust), clust_e50(i_pc, :, j_clust));
            interepoch50corr_4(i_clust, j_clust, i_pc) = tmp(1,2);
        end
    end
end



figure
tiledlayout(5,3);    
colors = ["red", "black", "green", "blue"];
for roi = 1:15
    nexttile;    
    imagesc(intercorr(:, :, roi));                        
    title('ROI PCA'+string(roi))
    colormap(gca,'parula');
    colorbar();
    caxis([-1 1])    
end




figure
tiledlayout(5,3);    
colors = ["red", "black", "green", "blue"];
for roi = 1:15
    nexttile;    
    imagesc(interepochcorr_second(:, :, roi));                        
    title('ROI PCA'+string(roi))
    colormap(gca,'parula');
    colorbar();
    caxis([-1 1])    
end



figure
tiledlayout(5,3);    
colors = ["red", "black", "green", "blue"];
for roi = 1:15
    nexttile;    
    imagesc(interepoch50corr_4(:, :, roi));                        
    title('ROI PCA'+string(roi))
    colormap(gca,'parula');
    colorbar();
    caxis([-1 1])    
end

save("FilteredBold_andEpoched");




%%

wind = 120;
epochs=1;
for nsub=indexsub
    read = 0;
    total = size(TS_W{1, nsub}, 2);
    i=0;
    while(total-read>wind)
        Wepochs(:, :, epochs) = TS_W{1, nsub}(:, i*wind+1:(i+1)*wind);
        i=i+1;
        read = i*wind;
        epochs=epochs+1;
    end
end

epochs=1;
for nsub=indexsub
    read = 0;
    total = size(TS_N3{1, nsub}, 2);
    i=0;
    while(total-read>wind)
        N3epochs(:, :, epochs) = TS_N3{1, nsub}(:, i*wind+1:(i+1)*wind);
        i=i+1;
        read = i*wind;
        epochs=epochs+1;
    end
end

for ep=1:size(Wepochs, 3)
    WFCepochs(:, :, ep) = corrcoef(squeeze(Wepochs(:, :, ep))');
end

for ep=1:size(N3epochs, 3)
    N3FCepochs(:, :, ep) = corrcoef(squeeze(N3epochs(:, :, ep))');
end
%%
total_success=0;
total_fail=0;
success_rois = [];
for i=1:90
    discarded=0;
    for j=1:90
        if j<=i
            continue
        end        
        %[p, ~, stats] = anova1([squeeze(WFCepochs(i, j, :)) squeeze(N3FCepochs(i,j,1:40))], [],'off');
        p = ranksum( squeeze(WFCepochs(i, j, :)), squeeze(N3FCepochs(i,j,1:40)));
        %tbl = multcompare(stats,"CriticalValueType","bonferroni");
        %if tbl(end)>0.05    
        %    discarded = discarded + 1;
        %end
        if p>(0.05/4005)
            discarded = discarded + 1;
        end
        
    end
    if discarded>0
        disp("Roi" + i+  " discarded " + discarded+ " correlations")
        total_fail = total_fail + 1;
    else
        disp("Roi" + i+  " did not fail")
        total_success = total_success + 1;
        success_rois(end+1) = i;
    end
end






