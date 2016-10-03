function [hp] = visgmm(obj);
% function [hp] = visgmm(obj);
mu = double(obj.mu);
sig = double(squeeze(obj.Sigma)');
hold on

for i = 1:length(mu);
    [x2,y2,z2] = ellipsoid(mu(i,1),mu(i,2),mu(i,3),sig(i,1),sig(i,2),sig(i,3),20);
    facecol = double(rand([1,3]));
    hf = surf(x2,y2,z2,'Edgecolor',facecol,'FaceColor',facecol,'EdgeAlpha',1, ...
              'FaceAlpha',0.5);
end
axis vis3d
rotate3d on
