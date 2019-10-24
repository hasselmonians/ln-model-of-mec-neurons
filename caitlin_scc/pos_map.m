function [posgrid, bins] = pos_map(pos, nBins, boxSize)

% take the histogram

boxUnit = boxSize / nBins;
bins = boxUnit/2:boxUnit:boxSize-boxUnit/2;

% store grid
posgrid = zeros(length(pos), nBins^2);

% loop over positions
for idx = 1:length(pos)

    % figure out the position index
    [~, xcoor] = min(abs(pos(idx,1)-bins));
    [~, ycoor] = min(abs(pos(idx,2)-bins));

    bin_idx = sub2ind([nBins, nBins], nBins - ycoor + 1, xcoor);

    posgrid(idx, bin_idx) = 1;

end

end
