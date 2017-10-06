function [cmap, kL, kU] = zerocmap(cmap)
% Makes the "zero" index of the current colormap white. (This need not be index 1.)

[cmin, cmax] = caxis();
L = length(cmap);

k = ceil(-cmin/(cmax-cmin)*L); % the "zero" index
if k<1, k=1; elseif k>L-1, k=L-1; end % clip
cmap(k,:) = 1;

% what's the range of data values assigned to the "zero" index k?
kL = (k-1)/L*(cmax-cmin)+cmin;
kU = k/L*(cmax-cmin)+cmin;

end