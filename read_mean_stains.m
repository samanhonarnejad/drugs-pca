function x = read_mean_stains(foldername, row, col, channel)
x = read_stains(foldername, row, col, channel);
if isempty(x)
    fprintf('warning: failed to read %d, %d from %s\n', row, col, foldername);
    return
end
x = mean(x);
