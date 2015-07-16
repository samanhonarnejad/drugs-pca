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
drug_effects = zeros(n_drugs, n_doses - 1, idx_pRB);
ch_647 = 3;
ch_568 = 5;

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
        folder = sprintf('%s_%s_%s_P%d/', cell_line, ...
            drugs{drug_id}, timepoint, plate_id);
        % DMSO was added to row 1.
        mn_ctrl = read_stain_reps(folder, 1, [col, col + 1], ch_647);
        for conc = 1 : n_doses - 1
            mn_drug = read_stain_reps(folder, 1 + conc, [col, col + 1], ...
                ch_647); 
            drug_effects(drug_id, conc, stain_id) = mn_drug ./ mn_ctrl;
        end
        col = col + 2;
    end
    % augment effect matrix by drug effects on p-Rb.
    mn_drug = zeros(n_doses, 1);
    for plate_id = 1 : 3
        folder = sprintf('%s_%s_%s_P%d/', cell_line, drugs{drug_id}, ...
            timepoint, plate_id);
        for conc = 1 : n_doses
            mn_drug(conc) = mn_drug(conc) + read_stain_reps(folder, ...
                conc, 1 : 12, ch_568);
        end
    end
    for conc = 2 : n_doses
        drug_effects(drug_id, conc - 1, idx_pRB) = mn_drug(conc) ./ ...
            mn_drug(1);
    end
    fprintf('finished analyzing stains for %s\n', drugs{drug_id});
end
save 10a_24h_drug_effects.mat

%%
n_sig = length(signals_647);
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
        folder = sprintf('%s_%s_%s_P%d/', cell_line, ...
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

%% PCA
load 10a_24h_drug_effects.mat
drug_effects = log2(drug_effects);
for k = 1 : n_sig + 1
    drug_effects(:, :, k) = drug_effects(:, :, k) ./ ...
        repmat(sd_stain(k), n_drugs, n_doses - 1);
end
n_drugs = length(drugs);
n_sigs = length(signals_647) + 1;
dfx = zeros(n_drugs * (n_doses - 1), n_sigs);
for drug_id = 1 : n_drugs
    dfx((1 : n_doses - 1) + (drug_id - 1) * (n_doses - 1), :) = ...
        squeeze(drug_effects(drug_id, :, :));
end
mn_dfx = mean(dfx);
dfx = dfx - repmat(mn_dfx, n_drugs * (n_doses - 1), 1);
coeff = pca(dfx);
offset = -mn_dfx * coeff;
sty = {'-ob', '-ok', '-om', '-^m', '-sm', '-or', '-^b', '-^k', '-^r'};
for drug_id = 1 : length(drugs)
    titr = squeeze(drug_effects(drug_id, :, :)) * coeff;
    x = titr(:, 1) - offset(1);
    y = titr(:, 2) - offset(2);
    plot(x, y, sty{drug_id}), hold('on');
end
idx_stain = find(sqrt(sum(coeff(:, [1, 2]) .^ 2, 2)) > 0.4);
for sig = idx_stain'
    vec = zeros(19, 1);
    vec(sig) = 0.05;
    vec = vec' * coeff;
    h = annotation('arrow');
    set(h, 'parent', gca(), 'position', [0, 0, vec(1), vec(2)]);
    text(1.2 * vec(1), 1.2 * vec(2), signals_647{sig});
end
xlim([-0.025, 0.075]);
ylim([-0.05, 0.05]);
legend(drugs);
hold('off');

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

%%
mn_cycif = mean(p27_cycif, 2);
mn_trad = mean(p27_trad, 2);
for d = 1 : 4
    
    line([p27_trad(d, 1), p27_trad(d, 2)], [mn_cycif(d), mn_cycif(d)]);
    line([mn_trad(d), mn_trad(d)], [p27_cycif(d, 1), p27_cycif(d, 2)]);
    
end
