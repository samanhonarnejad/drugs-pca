%% Initialize workspace.
drugs = { 'AZD6244', 'BEZ235', 'Erlotinib', 'Gefitinib', 'Lapatinib', ...
    'MK2206', 'PD0325901', 'PP242', 'Triciribine' };
% 18 signals from Alexa647 channels
signals_647 = { 'pRSK_{Ser380}', 'pERK_{Thr202,Tyr204}', 'pAKT_{Ser473}', ...
    'Foxo3a', 'pS6_{Ser235,Ser236}', 'p4EBP1_{Thr37,46}', 'pCDK2_{Tyr15}', ...
    'pCDK1_{Tyr15}', 'pP57_{Thr310}', 'pP27_{Ser10}', 'pP27_{Thr187}', ...
    'Survivin', 'p27', 'p21', 'p57', 'FoxM1', 'CyclinB', 'CyclinA' };
doses = [0, (10 * (3.1623 .^ (-6 : 0))) * 10 ^ -6]; % drug conc in M
n_doses = length(doses);
cell_line = '10a';
timepoint = '24h';
idx_pRB = length(signals_647) + 1;
n_drugs = length(drugs);
% drug_effects = nan(n_drugs, n_doses - 1, 2, idx_pRB);
drug_effects = struct('rep', {});
ch_647 = 3;
ch_568 = 5;
data_path = ['/Volumes/imstor/sorger/data/Operetta/Saman/Fixed Cell/', ...
    'Nature Scientific Data/Data/Drug Response/Analysis/'];
if ~exist(data_path, 'dir')
    error('Columbus export directory not mapped.');
end

%% Read mean of stainings for each drug treatment.
for drug_id = 1 : n_drugs
    % plate number starts over for each drug.
    plate_id = 0;
    for stain_id = 1 : length(signals_647)
        if mod(stain_id, 6) == 1
            % switch plates.
            plate_id = plate_id + 1;
            col = 1;
        end
        folder = sprintf([data_path, '%s_%s_%s_P%d/'], cell_line, ...
            drugs{drug_id}, timepoint, plate_id);
        % DMSO was added to row 1.
        mn_ctrl = mean(read_stain_reps(folder, 1, [col, col + 1], ch_647));
        for conc = 1 : n_doses - 1
            drug_reps = read_stain_reps(folder, 1 + conc, ...
                [col, col + 1], ch_647);
            n_rep = length(drug_reps);
            drug_effects(drug_id, conc, stain_id).rep = zeros(n_rep, 1);
            for rep = 1 : n_rep
                drug_effects(drug_id, conc, stain_id).rep(rep) = ...
                    drug_reps(rep) ./ mn_ctrl;
            end
        end
        col = col + 2;
    end
    % augment effect matrix by drug effects on p-Rb.
    for k = 1 : n_doses
        drug_effects(drug_id, conc, idx_pRB).rep = [];
    end
    for conc = 1 : n_doses
        mn_drug = [];
        for plate_id = 1 : 3
            folder = sprintf([data_path, '%s_%s_%s_P%d/'], cell_line, ...
                drugs{drug_id}, timepoint, plate_id);
            mn_drug = [mn_drug; read_stain_reps(folder, ...
                conc, 1 : 12, ch_568)]; %#ok<AGROW>
            n_rep = n_rep + 1;
        end
        if conc == 1
            mn_ctrl = mean(mn_drug(~isnan(mn_drug)));
        else
            drug_effects(drug_id, conc - 1, idx_pRB).rep = ...
                mn_drug(~isnan(mn_drug)) / mn_ctrl;
        end
    end
    fprintf('finished analyzing stains for %s\n', drugs{drug_id});
end
save 10a_24h_drug_effects.mat

%% PCA
load 10a_24h_drug_effects.mat
norm_dyn = false; %#ok<*UNRCH>
n_sig = length(signals_647);
if norm_dyn
    sd_stain = zeros(n_sig + 1, 1);
    for drug_id = 1 : n_drugs
        % plate number starts over for each drug.
        plate_id = 0;
        for stain_id = 1 : n_sig
            if mod(stain_id, 6) == 1
                % switch plates.
                plate_id = plate_id + 1;
                col = 1;
            end
            folder = sprintf([data_path, '%s_%s_%s_P%d/'], cell_line, ...
                drugs{drug_id}, timepoint, plate_id);
            % DMSO was added to row 1.
            sd_stain(stain_id) = sd_stain(stain_id) + ...
                std(read_stains(folder, 1, col, ch_647));
            col = col + 2;
        end
        sd_stain(idx_pRB) = sd_stain(idx_pRB) + ...
            std(read_stains(folder, 1, 6, ch_568));
    end
    sd_stain = sd_stain / n_drugs;
    arr_zero = [0.03, -0.03];
    arr_len = 0.025;
else
    sd_stain = ones(n_sig + 1, 1);
    arr_zero = [-1.4, -0.5];
    arr_len = 0.25;
end
mn_drug_effects = zeros(9, 7, 19);
for drug_id = 1 : 9
    for conc = 1 : 7
        for stain = 1 : 19
            drug_effects(drug_id, conc, stain).rep = log2(drug_effects( ...
                drug_id, conc, stain).rep) / sd_stain(stain);
            mn_drug_effects(drug_id, conc, stain) = ...
                mean(drug_effects(drug_id, conc, stain).rep);
        end
    end
end
figure(1), clf();
axes('position', [0.3, 0.3, 0.65, 0.65]);
n_drugs = length(drugs);
n_sigs = length(signals_647) + 1;
dfx = zeros(n_drugs * (n_doses - 1), n_sigs);
for drug_id = 1 : n_drugs
    dfx((1 : n_doses - 1) + (drug_id - 1) * (n_doses - 1), :) = ...
        squeeze(mn_drug_effects(drug_id, :, :));
end
mn_dfx = mean(dfx);
dfx = dfx - repmat(mn_dfx, n_drugs * (n_doses - 1), 1);
coeff = pca(dfx);
offset = -mn_dfx * coeff;
sty = {'-ob', '-ok', '-om', '-^m', '-sm', '-or', '-^b', '-^k', '-^r'};
shift_mat = repmat(mn_dfx, 7, 1);
h_curve = zeros(length(drugs), 1);
for drug_id = 1 : length(drugs)
    titr = (squeeze(mn_drug_effects(drug_id, :, :)) - shift_mat) * coeff;
    x = titr(:, 1) - offset(1);
    y = titr(:, 2) - offset(2);
    for conc = 1 : 7
        sd = zeros(1, 19);
        for stain = 1 : 19
            sd(stain) = std(drug_effects(drug_id, conc, stain).rep);
        end
        eb_mid_base = squeeze(mn_drug_effects(drug_id, conc, :))' - mn_dfx;
        eb_hi = (eb_mid_base + sd) * coeff - offset;
        eb_lo = (eb_mid_base - sd) * coeff - offset;
        line([eb_hi(1), eb_lo(1)], [y(conc), y(conc)]);
        line([x(conc), x(conc)], [eb_hi(2), eb_lo(2)]);
    end
    h_curve(drug_id) = plot(x, y, sty{drug_id});
    hold('on');
end
pca_compass(coeff, signals_647, arr_zero, arr_len, 0.4);
legend(h_curve, drugs);
hold('off');
% y-axis, second component
axes('position', [0.05, 0.3, 0.2, 0.65]);
bar(coeff(:, 2));
set(gca(), 'xtick', 1 : n_sig + 1, 'xticklabel', [signals_647, {'p-Rb'}]);
view(-90, 90);
% x-axis, first component
axes('position', [0.3, 0.05, 0.65, 0.2]);
bar(coeff(:, 1));
set(gca(), 'xtick', 1 : n_sig + 1, 'xticklabel', [signals_647, {'p-Rb'}]);

%% Correlation of drug effects on p27 with CycIF dataset
% traditional dataset
p27_trad = zeros(n_drugs, 2);
idx_p27 = find(strcmp(signals_647, 'p27'));
conc_corr = 6;
for drug_id = 1 : n_drugs
    % plate number starts over for each drug.
    plate_id = 0;
    for stain_id = 1 : length(signals_647)
        if mod(stain_id, 6) == 1
            % switch plates.
            plate_id = plate_id + 1;
            col = 1;
        end
        if stain_id == idx_p27
            folder = sprintf('%s_%s_%s_P%d/', cell_line, ...
                drugs{drug_id}, timepoint, plate_id);
            mn_ctrl = read_mean_stains(folder, 1, col, ch_647);
            mn_drug = read_mean_stains(folder, 1 + conc_corr, ...
                col, ch_647); 
            p27_trad(drug_id, 1) = log2(mn_drug ./ mn_ctrl);
            mn_ctrl = read_mean_stains(folder, 1, col + 1, ch_647);
            mn_drug = read_mean_stains(folder, 1 + conc_corr, ...
                col + 1, ch_647);
            p27_trad(drug_id, 2) = log2(mn_drug ./ mn_ctrl);
        end
        col = col + 2;
    end
end
idx_selu = find(strcmp(drugs, 'AZD6244'));
idx_bez = find(strcmp(drugs, 'BEZ235'));
idx_lapa = find(strcmp(drugs, 'Lapatinib'));
idx_pp242 = find(strcmp(drugs, 'PP242'));
p27_trad = p27_trad([idx_lapa, idx_selu, idx_bez, idx_pp242], :);
% cycif dataset
p27_cycif = zeros(4, 2);
for col = first_col : last_col
    if col == first_col
        ctrl = loadcycif(2, col, 'exclude', ignore);
    else
        other_ctrl = loadcycif(2, col, 'exclude', ignore);
        ctrl.data = [ctrl.data; other_ctrl.data];
    end
end
mn_ctrl = mean(ctrl.data(:, idx_p27));
idx_p27 = find(strcmp(ctrl.names, 'p27'));
conc_corr = 5;
for d = 1 : 4
    for r = 1 : 2
        rep = loadcycif(conc_corr, drug(d).col(r), 'exclude', ignore);
        mn_drug = mean(rep.data(:, idx_p27));
        p27_cycif(d, r) = log2(mn_drug ./ mn_ctrl);
    end
end
mn_cycif = mean(p27_cycif, 2);
mn_trad = mean(p27_trad, 2);
for d = 1 : 4
    line([p27_trad(d, 1), p27_trad(d, 2)], [mn_cycif(d), mn_cycif(d)]);
    line([mn_trad(d), mn_trad(d)], [p27_cycif(d, 1), p27_cycif(d, 2)]);
end

%% Histograms of pCDK2 and p27 in response to gefitinib.
n_bin = 100;
idx_pcdk2 = find(strcmp(signals_647, 'pCDK2_{Tyr15}'));
plate_id = ceil(idx_pcdk2 / 6);
folder = sprintf('%s_%s_%s_P%d/', cell_line, 'Gefitinib', timepoint, ...
    plate_id);

col1 = 2 * mod(7 - 1, 6) + 1;
col2 = col1;

pcdk2 = cell(1, n_doses);
pcdk2_lo = inf;
pcdk2_hi = -inf;
for conc = 1 : n_doses
    rep1 = read_stains(folder, conc, col1, ch_647);
    rep2 = read_stains(folder, conc, col2, ch_647);
    pcdk2{conc} = [rep1; rep2];
    % identify range
    min_pcdk2 = min(pcdk2{conc});
    if min_pcdk2 < pcdk2_lo
        pcdk2_lo = min_pcdk2;
    end
    max_pcdk2 = prctile(pcdk2{conc}, 99);
    if max_pcdk2 > pcdk2_hi
        pcdk2_hi = max_pcdk2;
    end
end

log_rng_cdk2 = linspace(log2(pcdk2_lo), log2(pcdk2_hi), n_bin);
figure();
for conc = 1 : n_doses
    n = hist(log2(pcdk2{conc}), log_rng_cdk2);
    n = n / sum(n);
    plot(log_rng_cdk2(1 : end - 1), n(1 : end - 1), 'k');
    hold on;
end
hold off;
