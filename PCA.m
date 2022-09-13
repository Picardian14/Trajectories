%% Preprocessing
%tsdata = zscore(tsdata, 0, 3);
%Wdata = zscore(Wdata, 0, 3);
%N3data = zscore(N3data, 0, 3);

[coeff,score,latent,tsq, explained, mu] = pca(mean(tsdataF,3)');
[coeff_roi,score_roi,latent_ror,tsq_roi, explained_roi, mu_roi] = pca(squeeze(mean(tsdataF, 1)));
figure
imagesc(mean(tsdata, 3))
title('original data')

figure
imagesc(mean(tsdataF, 3))
title('original data')
%%
% tsdata_pca has the subjects in the last dimension. This means that it
% will cluster subjects acording to a space of ROIs which explain most of
% the brain data
NCOMP = [1:15]; % 93% variance explained
clear tsdata_pca;
for i=1:size(tsdataF, 3)    
    tsdata_pca(:, :, i) = coeff(:, NCOMP)'*tsdataF(:, :, i);
end

% tsdata_pca_roi has the ROIS in the last dimension. This means that it
% will cluster ROIs in a reduced by PCA subject space
NCOMP = [1:5]; %92% variance explained
for i=1:size(tsdata, 1)
    tsdata_pca_roi(i, :, :) = coeff_roi(:, NCOMP)'*squeeze(tsdata(i, :, :))';
end
tsdata_pca_roi = permute(tsdata_pca_roi, [2 3 1]);
figure
imagesc(mean(tsdata_pca, 3))
title('PCA data')


%%
%{
figure
imagesc(mean(tsdata_reconstructed, 3))
title('Data reconstructed from PCA')

figure
imagesc(explained)
title('Variance explained')
%}
[Wcoeff,Wscore,Wlatent,Wtsq, Wexplained, Wmu] = pca(mean(Wdata,3)');
[N3coeff,N3score,N3latent,N3tsq, N3explained, N3mu] = pca(mean(N3data,3)');

for i=1:size(Wdata, 3)
    Wdata_pca(:, :, i) = coeff(:, NCOMP)'*Wdata(:, :, i);
end

for i=1:size(N3data, 3)
    N3data_pca(:, :, i) = coeff(:, NCOMP)'*N3data(:, :, i);
end
