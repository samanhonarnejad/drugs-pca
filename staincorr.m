function [rho_y, names] = staincorr(rows, cols, idx_y)

rho_y = [];
for row = rows(1) : rows(2)
  for col = cols(1) : cols(2)
    well = loadcycif(2, col, 'normalize', true, 'exclude', ...
        [16 : 24, 33, 34]);
    n_stains = size(well.data, 2) - 1;
    if isempty(rho_y)       
        rho_y = zeros(rows(2) - rows(1) + 1, cols(2) - cols(1) + 1, ...
            n_stains);
    end
    % measure correlation with all other stains
    rho = corr(well.data);
    rho_y_others = rho(idx_y, :);
    rho_y_others(idx_y) = [];
    rho_y(row-rows(1) + 1, col - cols(1) + 1, :) = rho_y_others;
  end
end
names = well.names;
names(idx_y) = [];
