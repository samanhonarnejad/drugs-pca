function well = loadcycif(row, col, varargin)
persistent pack names;
if isempty(pack)
    pack = load('matlab-CycIF-Saman-20150522', 'wellsum');
end
bad = [1, 8, 9, 32];
if isempty(names)
    [~, ~, raw] = xlsread('Saman-CycIF-metadata-20150521.xlsx', 'Channel');
    names = raw(:, 2);
    names(bad) = [];
end
exclude = [];
for k = 1 : 2: length(varargin)
    switch varargin{k}
        case 'exclude'
            exclude = varargin{k + 1};
        otherwise
            error(['unknown argument: ', varargin{k}]);
    end
end
data = pack.wellsum{row, col};
data(:, bad) = [];
data(:, exclude) = [];
well.data = data;
well.names = names;
well.names(exclude) = [];
