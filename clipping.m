function [yrng, xrng] = clipping(data, yrng0, xrng0)
xrng = prctile(data(:, 2), [0.01, 99.99]);
xctr = mean(xrng);
xrng = (xrng - xctr) * 1.1 + xctr;
xrng = linspace(xrng(1), xrng(2), 100);
yrng = prctile(data(:, 1), [0.01, 99.99]);
yctr = mean(yrng);
yrng = (yrng - yctr) * 1.1 + yctr;
yrng = linspace(yrng(1), yrng(2), 100);
if nargin == 3
    if xrng0(1) < xrng(1)
        xrng = linspace(xrng0(1), xrng(end), 100);
    end
    if xrng0(2) > xrng(end)
        xrng = linspace(xrng(1), xrng0(2), 100);
    end
    if yrng0(1) < yrng(1)
        yrng = linspace(yrng0(1), yrng(end), 100);
    end
    if yrng0(2) > yrng(end)
        yrng = linspace(yrng(1), yrng0(2), 100);
    end    
end
