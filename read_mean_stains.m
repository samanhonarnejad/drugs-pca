function x = read_mean_stains(foldername, row, col, channel)
x = mean(read_stains(foldername, row, col, channel));
