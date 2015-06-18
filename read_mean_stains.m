function x = read_mean_stains(foldername, row, col, channel)
if foldername(end) ~= '/'
    foldername = [foldername, '/'];
end
wellname = ['A' - 1 + row, num2str(col)];
csv_pattern = ['*result.', wellname, '[*Selected[1].csv'];
L = dir([foldername, csv_pattern]);
if isempty(L)
    error(['No CSV file in folder ', foldername, ' for well ', wellname]);
end
filename = [foldername, L.name];
csv_content = dlmread(filename, ',', 1, 16);
x = mean(csv_content(:, channel));
