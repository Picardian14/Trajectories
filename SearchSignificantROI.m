function success_rois = SearchSignificantROI(TS_W,TS_N3, indexsub)

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


