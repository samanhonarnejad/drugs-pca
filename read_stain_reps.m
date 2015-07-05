function mn_sig = read_stain_reps(folder, row, col, channel)
mn_sig = 0;
n = 0;
parfor rep = col
    mn_stain = read_mean_stains(folder, row, rep, channel);
    if ~isempty(mn_stain)
        mn_sig = mn_sig + mn_stain;
        n = n + 1;
    end
end
if n > 0
    mn_sig = mn_sig / n;
else
    error('No CSV file for any replicate in %s for row %d.', folder, row);
end
