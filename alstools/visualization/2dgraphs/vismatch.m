function [] = vismatch(trees,chm);
% function [] = vismatch(trees,chm);
% visualizes matched trees on CHM

siz = 8;
ii = ~isnan(trees.alsx);
imagesc(chm.x,chm.y,chm.z);hold on
plot(trees.x(ii),trees.y(ii),'.g','markersize',siz);
plot(trees.alsx(ii),trees.alsy(ii),'.r','markersize',siz);
hl = line([trees.alsx,trees.x]',[trees.alsy,trees.y]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
axis([min(chm.x) max(chm.x) min(chm.y) max(chm.y)]);
axis equal;axis tight;axis xy;
colormap(gray)
hcb = colorbar('East');
axes(hcb);ylabel('Canopy Height [m]');set(gca,'yaxislocation','left','color','w','ycolor','w','xcolor','w');
set(gcf,'inverthardcopy','off','color','w');
setstrings(12)