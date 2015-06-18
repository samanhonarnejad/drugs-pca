%% Initialize workspace.
drugs = { 'AZD6244', 'BEZ235', 'Erlotinib', 'Gefitinib', 'Lapatinib', ...
    'MK2206', 'PD0325901', 'PP242', 'Triciribine' };
% 18 signals from Alexa647 channels
signals_647 = { 'pRSK_{Ser380}', 'pERK_{Thr202,Tyr204}', 'pAKT_{Ser473}', ...
    'Foxo3a', 'pS6_{Ser235,Ser236}', 'p4EBP1_{Thr37,46}', 'pCDK2_{Tyr15}', ...
    'pCDK1_{Tyr15}', 'pP57_{Thr310}', 'pP27_{Ser10}', 'pP27_{Thr187}', ...
    'Survivin', 'p27', 'p21', 'p57', 'FoxM1', 'CyclinB', 'CyclinA' };
doses = [(10 * (3.1623.^(0:-1:-6)))*10^-6 0]; % drug concentration in M
doses = doses(8:-1:1);
cell_line = '10a';
timepoint = '24h';
idx_pRB = length(signals_647) + 1;
drug_effects = zeros(length(drugs), idx_pRB);
row_ctrl = 1;
row_drug = 8;
ch_647 = 3;
ch_568 = 5;

%% Read mean of stainings for each drug treatment.
for drug_id = 1 : length(drugs)
    % plate number starts over for each drug.
    plate_id = 0;
    for stain_id = 1 : length(signals_647)
        if mod(stain_id, 6) == 1
            % switch plates.
            plate_id = plate_id + 1;
            col = 0;
        end
        foldername = sprintf('%s_%s_%s_P%d/', cell_line, ...
            drugs{drug_id}, timepoint, plate_id);
        col = col + 1;
        % highest dose was applied to row 8, DMSO ro row 1. 
        mn_ctrl = read_mean_stains(foldername, row_ctrl, col, ch_647);
        mn_drug = read_mean_stains(foldername, row_drug, col, ch_647);
        col = col + 1;
        mn_ctrl = (mn_ctrl + read_mean_stains(foldername, row_ctrl, ...
            col, ch_647)) / 2;
        mn_drug = (mn_drug + read_mean_stains(foldername, row_drug, ...
            col, ch_647)) / 2;
        drug_effects(drug_id, stain_id) = mn_drug ./ mn_ctrl;
    end
    % augment effect matrix by drug effects on p-Rb.
    mn_drug = 0;
    mn_ctrl = 0;
    for plate_id = 1 : 3
        foldername = sprintf('%s_%s_%s_P%d/', cell_line, ...
            drugs{drug_id}, timepoint, plate_id);
        for col = 1 : 12
            mn_drug = mn_drug + read_mean_stains(foldername, ...
                row_drug, col, ch_568);
            mn_ctrl = mn_ctrl + read_mean_stains(foldername, ...
                row_ctrl, col, ch_568);
        end
    end
    drug_effects(drug_id, idx_pRB) = mn_drug ./ mn_ctrl;
end
save 10a_24h_drug_effects.mat

%% PCA
load 10a_24h_drug_effects.mat
n_drugs = length(drugs);
dfx = log2(drug_effects);
dfx = dfx - repmat(mean(dfx), n_drugs, 1);
[coeff, effects_pca, latent] = pca(dfx);
x = effects_pca(:,1);
y = effects_pca(:,2);
clf();
plot (x, y, '*');
for k = 1 : length(drugs)
    text(x(k), y(k), drugs{k});
end
idx_stain = [];
for c = 1 : 2
    idx_stain = unique([idx_stain; find(coeff(:, c) > 0.4)]);
end
for a = idx_stain'
    vec = zeros(19, 1);
    vec(a) = 1;
    vec = vec' * coeff;
    h = annotation('arrow');
    set(h, 'parent', gca(), 'position', [0, 0, vec(1), vec(2)]);
    text(1.1 * vec(1), 1.1 * vec(2), signals_647{a});
end
