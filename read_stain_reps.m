function sig = read_stain_reps(folder, row, col, channel)
sig = nan(length(col), 1);
n = 0;
for rep = col
    mn_stain = read_mean_stains(folder, row, rep, channel);
    if ~isempty(mn_stain)
        n = n + 1;
        sig(n) = mn_stain;
    end
end
if n == 0
    error('No CSV file for any replicate in %s for row %d, col %d, %d', ...
        folder, row, col(1), col(2));
end
sig = sig(1 : n);
