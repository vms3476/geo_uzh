function [mod] = subsetmod(mod,are);
% function [mod] = subsetmod(mod,are);
% subset model structure to region of interest defined by are

if isstruct(are)
    are = [min(are.x) max(are.x) min(are.y) max(are.y)];
end
ii = mod.x >= are(1) & mod.x <= are(2);
jj = mod.y >= are(3) & mod.y <= are(4);
mod.x = mod.x(ii);
mod.y = mod.y(jj);
mod.z = mod.z(jj,ii);
if isfield(mod,'X')
    mod.X = mod.X(jj,ii);
    mod.Y = mod.Y(jj,ii);
end