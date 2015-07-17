%%
load Organize_10a_Ligand_Response_Single_Cell_Data.mat;
idx_p27 = find(strcmp(signals_647, 'p27'));
idx_pS6 = find(strcmp(signals_647, 'pS6_Ser235_236'));

%% determine range for histograms
n_egf_doses = 8;
pRb_lo = inf;
pRb_hi = -inf;
p27_lo = inf;
p27_hi = -inf;
pS6_lo = inf;
pS6_hi = -inf;
for dose = 1 : n_egf_doses
    for r = 1 : 30
        pRb = squeeze(single_cell_Alexa568_nuc_data(1, 1, dose, r, 1, :));
        min_pRb = min(pRb);
        if min_pRb < pRb_lo
            pRb_lo = min_pRb;
        end
        max_pRb = prctile(pRb, 99);
        if max_pRb > pRb_hi
            pRb_hi = max_pRb;
        end
    end
    for r = idx_p27'
        p27 = squeeze(single_cell_Alexa647_nuc_data(1, 1, dose, r, 1, :));
        min_p27 = min(p27);
        if min_p27 < p27_lo
            p27_lo = min_p27;
        end
        max_p27 = prctile(p27, 99);
        if max_p27 > p27_hi
            p27_hi = max_p27;
        end
    end
    pS6 = squeeze(single_cell_Alexa647_nuc_data(1, 1, dose, idx_pS6, 1, :));
    min_pS6 = min(pS6);
    if min_pS6 < pS6_lo
        pS6_lo = min_pS6;
    end
    max_pS6 = prctile(pS6, 99);
    if max_pS6 > pS6_hi
        pS6_hi = max_pS6;
    end
end

%% calculate histograms
n_bin = 100;
n_rep = 30;
fn_log = @log2;
log_rng_pRb = linspace(fn_log(pRb_lo), fn_log(pRb_hi), n_bin);
log_rng_p27 = linspace(fn_log(p27_lo), fn_log(p27_hi), n_bin);
log_rng_pS6 = linspace(fn_log(pS6_lo), fn_log(pS6_hi), n_bin);
pRb = zeros(n_rep, n_bin, n_egf_doses);
p27 = zeros(length(idx_p27), n_bin, n_egf_doses);
pS6 = zeros(n_bin, n_egf_doses);
for dose = 1 : n_egf_doses
    for rep = 1 : n_rep
        rep_pRb = squeeze(single_cell_Alexa568_nuc_data(1, 1, dose, ...
            rep, 1, :));
        n = hist(fn_log(rep_pRb), log_rng_pRb);
        pRb(rep, :, dose) = n ./ sum(n);
    end
    for rep = 1 : length(idx_p27)
        rep_p27 = squeeze(single_cell_Alexa647_nuc_data(1, 1, dose, ...
            idx_p27(rep), 1, :));
        n = hist(fn_log(rep_p27), log_rng_p27);
        p27(rep, :, dose) = n' ./ sum(n);
    end
    rep_pS6 = squeeze(single_cell_Alexa647_nuc_data(1, 1, dose, ...
        idx_pS6, 1, :));
    n = hist(fn_log(rep_pS6), log_rng_pS6);
    pS6(:, dose) = n' ./ sum(n);
end
md_pRb = squeeze(median(pRb));
mn_p27 = squeeze(mean(p27));

%% plot histograms
n = 1;
for k = 2 : n_egf_doses
    subplot(7, 3, 3 * n - 2);
    plot(log_rng_pRb(1 : end - 1), md_pRb(1 : end - 1, 1), 'k');
    hold on;
    plot(log_rng_pRb(1 : end - 1), md_pRb(1 : end - 1, k), 'r');
    hold off;
    subplot(7, 3, 3 * n - 1);
    plot(log_rng_p27(1 : end - 1), mn_p27(1 : end - 1, 1), 'k');
    hold on;
    plot(log_rng_p27(1 : end - 1), mn_p27(1 : end - 1, k), 'r');
    hold off;
    subplot(7, 3, 3 * n);
    plot(log_rng_pS6(1 : end - 1), pS6(1 : end - 1, 1), 'k');
    hold on;
    plot(log_rng_pS6(1 : end - 1), pS6(1 : end - 1, k), 'r');
    hold off;
    n = n + 1;
end

%% 1D scatter plots
for k = 1 : n_egf_doses
    ss_pRb = fn_log(squeeze(single_cell_Alexa568_nuc_data(1, 1, k, 1, 1, :)));
    ss_pRb(isnan(ss_pRb)) = [];
    
    pRb_idx = ss_pRb - log_rng_pRb(1);
    pRb_idx = pRb_idx ./ (log_rng_pRb(end) - log_rng_pRb(1));
    pRb_idx = round(pRb_idx * (n_bin - 1) + 1);
    
    exc_pRb = pRb_idx > 100;
    ss_pRb(exc_pRb) = [];
    pRb_idx(exc_pRb) = [];

    xjit = randn(length(ss_pRb), 1) .* pRb(1, pRb_idx, k)';
    subplot(1, 8, k);
    plot(xjit, ss_pRb, '.');
    ylim([7, 11]);
end
