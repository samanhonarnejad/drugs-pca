%%
ignore = [5, 6, 14, 16 : 24, 30, 31, 33, 34];
drug = struct('name', {'Lap', 'Selu', 'Dact', 'PP242'}, 'col', {[4, 5], ...
    [6, 7], [8, 9], [10, 11]});

%% PCA
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
        delta_drug(n, :) = mn_drug ./ mn_ctrl - 1;
    end
end
% principal component analysis of drug-induced shifts
coeff = pca(delta_drug);
% plot projection along first two principal components
figure(1), clf(), hold('on');
col = {[0.3, 0.1, 1], [1, 0.1, 0.5], [0.6, 1, 0.0], [0.5, 0.2, 0]};
idx0 = 0 : 5 : 15;
for d = 1 : 4
    mn_prev = [];
    for c = 2 : 6
        rep1 = loadcycif(c, drug(d).col(1), 'exclude', ignore);
        rep1 = (mean(rep1.data) - mn_ctrl) ./ mn_ctrl * coeff;
        rep2 = loadcycif(c, drug(d).col(2), 'exclude', ignore);
        rep2 = (mean(rep2.data) - mn_ctrl) ./ mn_ctrl * coeff;
        mn = 0.5 * (rep1 + rep2);
        ellipse2d(rep1, rep2, col{d}, (c / 6) .^ 2);
        if ~isempty(mn_prev)
            line([mn_prev(1), mn(1)], [mn_prev(2), mn(2)], 'color', ...
                col{d}, 'linewidth', (c - 2) .^ 1.2);
        end
        mn_prev = mn;
    end
end

% 1: overall trend with p-RB
% 2: p27
% 3: p-S6

%%
% The PCA of titration effects show a dose for each drug that captures the
% most variation
conc_idx = [5, 3, 2, 3];
delta_drug = zeros(4, n_ch);
limit = zeros(8, n_ch);
for d = 1 : 4
    rep1 = loadcycif(conc_idx(d), drug(d).col(1), 'exclude', ignore);
    rep2 = loadcycif(conc_idx(d), drug(d).col(2), 'exclude', ignore);
    pooled = [rep1.data; rep2.data];
    limit((d * 2 - 1) : (d * 2), :) = prctile(pooled, [1, 99]);
    mn_drug = mean(pooled);
    n = n + 1;
    delta_drug(n, :) = mn_drug ./ mn_ctrl - 1;    
end
coeff = pca(delta_drug);
max_sig = max(limit);
max_sig_pca = max_sig * coeff;
figure();
for d = 1 : 4
    rep = loadcycif(conc_idx(d), drug(d).col(1), 'exclude', ignore);
    rep_pca = rep.data * coeff;
    subplot(1, 4, d);
    scatplot(rep_pca(:, 1), rep.data(:, 7));
    xlim([0, 1.1 * max_sig_pca(1)]);
    ylim([0, 1.1 * max_sig(7)]);
end

%% Cluster analysis of single-cell states
% identify clusters of signaling states in control well. Use cosine, i.e.
% angle between stain vectors as distance metric such that even one stain
% can make a difference. Begin with a large number of clusters k.
close all;
conc_idx = [3, 3, 4, 3];
drug_id = 1;
disp(drug(drug_id).name);
rep = loadcycif(conc_idx(drug_id), drug(drug_id).col(1), 'exclude', ignore);

% rep = loadcycif(2, 4, 'exclude', ignore);
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
        subplot(4, 4, p);
        p = p + 1;
        plot([x1; x2; x_both]', [n1; n2; n_both]');
        title(sprintf('cluster %d vs cluster %d', k1, k2));
    end
end

%% this was done by eye because cluster labels change
idx(idx == 2) = 1;
idx(idx == 5) = 2;
k = 4;
figure(); % ... go back to plotting histograms

%%
idx(idx == 3) = 1;
idx(idx == 4) = 3;
k = 3;
figure(); % ... go back to plotting histograms

%% Characterize clusters
% Eventually, we converge at three clusters, with one very small cluster of
% about 1.5% of the whole population. We calculate the vector between the
% center of that cluster and all other cells to see what makes the
% difference. No interesting observations here for any of the other drugs.

c1 = mean(rep.data(idx ~= 3, :));
c2 = mean(rep.data(idx == 3, :));
[shift_sorted, idx_stain] = sort(abs(c2 - c1), 'descend');

%%
disp(rep.names(idx_stain(1 : 3)));

% pS6(235) and pS6(240) distinguish a new cluster after Lapatinib
% treatment.
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




