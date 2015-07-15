%% Initialize workspace.
ignore = [5, 6, 14, 16 : 24, 30, 31, 33, 34];
drug = struct('name', {'Lap', 'Selu', 'Dact', 'PP242'}, 'col', {[4, 5], ...
    [6, 7], [8, 9], [10, 11]});

%% PCA of drug effects on bulk signaling states.
% load and average all controls
first_col = 4;
last_col = 11;
for col = first_col : last_col
    if col == first_col
        ctrl = loadcycif(2, col, 'exclude', ignore);
    else
        other_ctrl = loadcycif(2, col, 'exclude', ignore);
        ctrl.data = [ctrl.data; other_ctrl.data];
    end
end
mn_ctrl = mean(ctrl.data);
n_ch = length(mn_ctrl);
% load drug titrations and calculate population shifts
n = 0;
delta_drug = zeros(20, n_ch);
res = struct('rep', {});
for d = 1 : 4
    for c = 2 : 6
        rep1 = loadcycif(c, drug(d).col(1), 'exclude', ignore);
        rep2 = loadcycif(c, drug(d).col(2), 'exclude', ignore);
        mn_drug = mean([rep1.data; rep2.data]);
        n = n + 1;
        % delta = (drug - ctrl) / ctrl = drug / ctrl - 1
        delta_drug(n, :) = log2(mn_drug ./ mn_ctrl);
        % delta_drug(n, :) = mn_drug ./ mn_ctrl - 1;
    end
end
save cycif_drug_effects.mat

%% principal component analysis of drug-induced shifts

coeff = pca(delta_drug);
% plot projection along first two principal components
figure(1), clf(), hold('on');
col = {[0.3, 0.1, 1], [1, 0.1, 0.5], [0.6, 1, 0.0], [0.5, 0.2, 0]};
idx0 = 0 : 5 : 15;
for d = 1 : 4
    mn_prev = [];
    for c = 2 : 6
        rep1 = loadcycif(c, drug(d).col(1), 'exclude', ignore);
        rep1 = log2(mean(rep1.data) ./ mn_ctrl) * coeff;
        rep2 = loadcycif(c, drug(d).col(2), 'exclude', ignore);
        rep2 = log2(mean(rep2.data) ./ mn_ctrl) * coeff;
        mn = 0.5 * (rep1 + rep2);
        ellipse2d(rep1, rep2, col{d}, (c / 6) .^ 2);
        if ~isempty(mn_prev)
            line([mn_prev(1), mn(1)], [mn_prev(2), mn(2)], 'color', ...
                col{d}, 'linewidth', (c - 2) .^ 1.2);
        end
        mn_prev = mn;
    end
end
% 1: overall trend with p-RB, but the two strongest components are pS6
% 2: p27
% 3: again p-S6
% CDK2 was not part of the CyCIF stains but p27 and CDK2 can be assumed to
% be opposing signaling components (check for anticorration between p27 and
% CDK2 in non-CyCIF dataset).

%% Cluster analysis of single-cell states in control condition
% identify clusters of signaling states in control well. Use cosine, i.e.
% angle between stain vectors as distance metric such that even one stain
% can make a difference. Begin with a large number of clusters k.
rep = loadcycif(2, 4, 'exclude', ignore);
n_cells = size(rep.data, 1);
sd_rep = sqrt(var(rep.data));
rep.data = rep.data ./ repmat(sd_rep, n_cells, 1);
k = 5;
idx = kmeans(rep.data, k, 'Replicates', 10, 'Distance', 'cosine');

%% histogram along projection onto axis between cluster centers
% Use histograms to check pairwise against oversegmentation: if the sum of
% the two histograms is not bimodal, the cluster boundary is not justified.
p = 1;
for k1 = 1 : k - 1
    for k2 = (k1 + 1) : k
        c1 = mean(rep.data(idx == k1, :));
        c2 = mean(rep.data(idx == k2, :));
        sep = (c2 - c1)' / norm(c2 - c1);
        subpop1 = rep.data(idx == k1, :) * sep;
        subpop2 = rep.data(idx == k2, :) * sep;
        [n1, x1] = hist(subpop1, 100);
        disp(sum(n1));
        [n2, x2] = hist(subpop2, 100);
        disp(sum(n2));
        [n_both, x_both] = hist(rep.data(idx == k1 | idx == k2, :) * sep, 100);
        subplot(4, 3, p);
        p = p + 1;
        plot([x1; x2; x_both]', [n1; n2; n_both]');
        title(sprintf('cluster %d vs cluster %d', k1, k2));
    end
end

%% this was done by eye because cluster labels change
idx(idx == 4) = 1;
idx(idx == 5) = 4;
k = 4;
clf(); % ... go back to plotting histograms

%%
idx(idx == 2) = 1;
idx(idx == 4) = 2;
k = 3;
clf(); % ... go back to plotting histograms

%% Examine which signaling states distinguish each cluster.
load kmeans-merged.mat
whole = mean(rep.data);
shift = zeros(k, n_ch);
for n = 1 : k
    clst = mean(rep.data(idx == n, :));
    shift(n, :) = clst - whole;
    [~, idx_stain] = sort(abs(shift(n, :)), 'descend');
    fprintf('\ncluster %d: %d\n', n, sum(idx == n));
    for p = 1 : 3
        fprintf('%s: %f\n', rep.names{idx_stain(p)}, ...
            shift(n, idx_stain(p)));
    end
end

% cluster 1: 11476     cluster 2: 202      cluster 3: 2859
% p-RB: -0.460307      p-H3: 7.771845      p-RB: 1.686127
% pS6(240): -0.333169  p-Aurora: 7.421638  PCNA: 1.201503
% pS6(235): -0.305965  gH2ax: 3.796556     pS6(240): 1.193047
%
% We believe that cluster 2 corresponds to mitotic cells. The two other
% clusters are much more populated and correspond to a pRb/pS6-high or
% pRb/pS6-low state. It is difficult to find more than two clusters in
% any of the drug-treated cells.

%% Visualize unperturbed single-cell clusters in PCA of cluster centers.
coeff_clust = pca(shift);
clust_pca = rep.data * coeff_clust;
[yrng, xrng] = clipping(clust_pca);
clf();
sty = {'k', 'b', 'r'};
for n = 1 : k
    im_ctrl = double(hist3(clust_pca(n == idx, [1, 2]), {yrng, xrng}));
    im_ctrl = im_ctrl ./ size(clust_pca, 1);
    im_ctrl = conv2(im_ctrl, fspecial('gaussian', [3, 3], 0.4), 'same');
    contour(im_ctrl, 0.19 .^ (1 : 6), sty{n}), hold('on');
end
ticks = 10 : 10 : 100;
set(gca(), 'xtick', ticks, 'xticklabel', round(xrng(ticks), 1), ...
    'ytick', ticks, 'yticklabel', round(yrng(ticks), 1));

%% Visualize cluster changes due to drugs in PCA of bulk drug effects.
ctrl = loadcycif(2, 5, 'exclude', ignore);
n_cells = size(ctrl.data, 1);
ctrl.data = ctrl.data ./ repmat(sd_rep, n_cells, 1);
ctrl_pca = ctrl.data * coeff;
[yrng, xrng] = clipping(clust_pca);
for drug_id = 1 : 4
    for dose = 3 : 7
        rep = loadcycif(dose, drug(drug_id).col(1), 'exclude', ignore);
        n_cells = size(rep.data, 1);
        rep.data = rep.data ./ repmat(sd_rep, n_cells, 1);
        drug_pca = rep.data * coeff;
        [yrng, xrng] = clipping(clust_pca, yrng, xrng);
    end
end
im_ctrl = double(hist3(ctrl_pca(:, [1, 2]), {yrng, xrng}));
im_ctrl = im_ctrl ./ sum(im_ctrl(:));
im_drug = cell(4, 1);
for drug_id = 1 : 4
    figure(drug_id);
    for dose = 2 : 7
        rep = loadcycif(dose, drug(drug_id).col(1), 'exclude', ignore);
        n_cells = size(rep.data, 1);
        rep.data = rep.data ./ repmat(sd_rep, n_cells, 1);
        drug_pca = rep.data * coeff;
        im_drug = double(hist3(drug_pca(:, [1, 2]), {yrng, xrng}));
        im_drug = im_drug ./ sum(im_drug(:));
        subplot(2, 3, dose - 1);
        contour(im_ctrl, linspace(-.015, .015, 12), 'k'), hold('on');
        contour(im_drug, linspace(-.015, .015, 12), 'r'), hold('off');
        title(drug(drug_id).name);
        xlim([10, 60]);
        ylim([5, 40]);
    end
end

%%
conc_idx = [6, 4, 3, 4];
drug_id = 4;
ctrl = loadcycif(conc_idx(drug_id), drug(drug_id).col(1), 'exclude', ignore);
idx = kmeans(ctrl.data, 2, 'Replicates', 10, 'Distance', 'cosine');
n_cells = size(ctrl.data, 1);
ctrl.data = ctrl.data ./ repmat(sd_rep, n_cells, 1);
indep = [5, 6, 10, 11, 12, 13, 14, 15];

%%
n = 0;
for k1 = 1 : length(indep) - 1
    for k2 = (k1 + 1) : length(indep)
        n = n + 1;
        subplot(4, 7, n);
        subset = (rand(n_cells, 1) < 0.1); 
        gscatter(log2(ctrl.data(subset, k1)), log2(ctrl.data(subset, ...
            k2)), idx(subset), 'kbr','ov^',[],'off');
        xlabel(ctrl.names{indep(k1)});
        ylabel(ctrl.names{indep(k2)});
    end
end
    
% p53 strikingly bimodal in Lapatinib case.
%
%

%% Correlations
[rho, names] = staincorr([2, 7], [4, 11], idx_y);
rho_ctrl = squeeze(rho(1, :, :));
mn_ctrl = mean(squeeze(rho_ctrl));
rho_lap = squeeze(rho(6, 1 : 2, :));
mn_lap = mean(rho_lap);
rho_selu = squeeze(rho(6, 3 : 4, :));
mn_selu = mean(rho_selu);
rho_dact = squeeze(rho(6, 5 : 6, :));
mn_dact = mean(rho_dact);
rho_PP242 = squeeze(rho(6, 5 : 6, :));
mn_PP242 = mean(rho_PP242);

figure(3), bar([mn_ctrl;mn_lap;mn_selu;mn_dact;mn_PP242]')
set(gca(), 'xtick', 1 : length(mn_ctrl), 'xticklabel', names);


%% compare correlations across drugs
n = 1;
for row = [2, 6]
    for col = 4 : 11
        subplot(2, 8, n);
        well = loadcycif(row, col, 'normalize', true, 'exclude', ...
            [16 : 24, 33, 34]);
        rho = corr(well.data);
        rho(eye(size(rho, 1)) == 1) = 0;
        imagesc(rho);
        n = n + 1;
    end
end


%%
% average control 
rho = [];
n = 0;
for col = 4 : 11
    well = loadcycif(2, col, 'normalize', true, 'exclude', ignore);
    if isempty(rho)
        rho = corr(well.data);
    else
        rho = rho + corr(well.data);
    end
    n = n + 1;
end
rho_ctrl = rho ./ n;

for n = 1 : 4
    drugrep1 = loadcycif(5, drug(n).col(1), 'normalize', true, 'exclude', ignore);
    drugrep2 = loadcycif(5, drug(n).col(2), 'normalize', true, 'exclude', ignore);
    rho_drug = (corr(drugrep1.data) + corr(drugrep2.data)) / 2;

    rho_diff = rho_drug - rho_ctrl;
    msk = tril(ones(size(rho, 1))) == 1;
    rho_diff(msk) = nan;

    mn_rho_diff = mean(rho_diff(~msk(:)));
    sd_rho_diff = std(rho_diff(~msk(:)));
    z_rho_diff = (rho_diff(:) - mn_rho_diff) / sd_rho_diff;
    subplot(1, 4, n);
    hist(z_rho_diff(~msk(:)), 50);
    title(drug(n).name);
    [z_rho_diff, idx] = sort(abs(z_rho_diff), 'descend');
    idx_z0 = find(~isnan(z_rho_diff), 1);

    fprintf('\n===== %s =====\n', drug(n).name);
    for k = idx_z0 : length(z_rho_diff)
        [row, col] = ind2sub(size(rho_diff), idx(k));
        z = z_rho_diff(k);
        fprintf('%0.1f\t%0.3f\t%0.3f\t%s\t%s\n', z, rho_ctrl(row, col), ...
            rho_diff(row, col), well.names{row}, well.names{col});
        if (abs(z) < 2)
            break
        end
    end

end


%%
figure(1);
y_pred = X_pca * beta;
subplot(1, 2, 1), plot(y_pred, y, '*');
title(num2str(corr(y_pred, y)));
% visualize 'importance' of different dimensions in X_pca
subplot(1, 2, 2), imagesc(beta);


%%
% we see that expcted individual readouts, like area and Ki-67 are good
% predictors.
[~, idx_dim] = max(abs(beta));
[~, idx_stain] = max(abs(coeff(:, idx_dim)));
disp(lbl_stains{idx_stain});

% conclusion: we can show proof of principle, but much easier to understand
% with simple correlation of measurements relative to p-Rb
% conclusion: since area and Ki-67 seem rather obvious, use PLSR to
% collapse predicted variables into a 'proliferation-readout' component.


%%
rho = corr(ctrl);
z = linkage(abs(rho));
n_clust = 25;
c = cluster(z, 'maxclust', n_clust);
n_clust = max(c);

rhoc = zeros(size(rho));
row = 1;
for dim = 1 : n_clust
    for src = find(c == dim)'
        rhoc(row, :) = rho(src, :);
        row = row + 1;
    end
end
rho = rhoc;
col = 1;
for dim = 1 : n_clust
    for src = find(c == dim)'
        rhoc(:, col) = rho(:, src);
        col = col + 1;
    end
end
figure(2), imagesc(rhoc);


%%

ctrl1 = loadcycif(2, 4, 'normalize', true);
idx_pRB = strmatch('p-RB', ctrl1.names);
rho1 = corr(ctrl1.data);
rho_pRB1 = rho1(:, idx_pRB);

ctrl2 = loadcycif(2, 5, 'normalize', true);
rho2 = corr(ctrl2.data);
rho_pRB2 = rho2(:, idx_pRB);

plot(rho_pRB1, rho_pRB2, '*');
% for k = 1 : length(rho_pRB1)
%     text(rho_pRB1(k), rho_pRB
% end
text(rho_pRB1+0.025, rho_pRB2-0.025, ctrl1.names);


%%

[rho_pRB, idx] = sort(rho_pRB);
figure(3);

for k = 1 : length(rho_pRB)
    line([rho_pRB(k), rho_pRB(k)], [-1, 1]);
    text(rho_pRB(k), 1.5, lbl_stains(idx(k)));
end
ylim([0, 2]);

% ugly! instead plot one replicate against the other




