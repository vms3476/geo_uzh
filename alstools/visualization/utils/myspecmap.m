function [cmap] = myspecmap(n);
%  function [cmap] = myspecmap(n);
% Colormap with Spectrum, builds upon createSpectrum.m
  
if nargin == 0
  n = 256;
end

[L,RGB] = createSpectrum('1964_FULL');

ocmap = squeeze(RGB(1,20:370,:));
nx = 1:n;
for i = 1:3
  cmap(:,i) = interp1(1:length(ocmap),ocmap(:,i),nx);
end
