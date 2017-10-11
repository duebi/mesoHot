function cmap = redbluecmap()
% Returns a red-blue colormap array. The cmap length is *not* 64.

cmap1 = colormap('hot');
cmap2 = colormap('jet');

red = cmap1(5:24,:); % dark to light
blue = cmap2(20:-1:1,:); % light to dark

% "stretch" out the cmap a bit
red2 = interp1(1:length(red),red,1:0.5:length(red));
blue2 = interp1(1:length(blue),blue,1:0.5:length(blue));

cmap = [red2; 1 1 1; blue2];


end