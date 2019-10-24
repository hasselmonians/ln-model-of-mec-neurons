function [posgrid, xBins, yBins] = pos_map(pos, nBins, boxSize)

% compute the size of a grid unit
boxUnit = boxSize / nBins;
% collect the center points of all the bins
xBins = boxUnit/2:boxUnit:xBoxSize-boxUnit/2;
yBins = boxUnit/2:boxUnit:yBoxSize-boxUnit/2;
% store grid
posgrid = zeros(length(pos), length(xBins) * length(yBins));

% compute the bin index for each data point
for idx = 1:length(pos)

    % figure out the position index
    [~, xcoor] = min(abs(pos(idx, 1) - xBins));
    [~, ycoor] = min(abs(pos(idx, 2) - yBins));

    bin_idx = sub2ind([nBins, nBins], nBins - ycoor + 1, xcoor);
    posgrid(idx, bin_idx) = 1;

end
