function x = read_mean_stains(foldername, row, col, channel)
x = read_stains(foldername, row, col, channel);
if isempty(x)
    return
end
x = mean(x);
