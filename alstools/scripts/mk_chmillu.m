if 0
load /Users/morsdorf/Desktop/ForestGrowth/change.mat
are = [812820    814790    170930    171840];
ii = ras.x >= are(1) & ras.x <= are(2);
jj = ras.y >= are(3) & ras.y <= are(4);

close all
figure(1);
shadmod(ras.x(ii),ras.y(jj),ras.dtm10(jj,ii));
swisstick;box on
colormap(terrain);
colorbar
print -depsc2 -r600 /Users/morsdorf/Graphiken/dtm_snp.eps

figure(2);
shadmod(ras.x(ii),ras.y(jj),ras.dsm10(jj,ii));
swisstick;box on
colormap(terrain);
colorbar
print -depsc2 -r600 /Users/morsdorf/Graphiken/dsm_snp.eps
end
figure(3);
myimage(ras.x(ii),ras.y(jj),ras.chm10(jj,ii));
swisstick;box on
colormap(vegetation2);
caxis([0 30])
colorbar
%print -depsc2 -r600 /Users/morsdorf/Graphiken/chm_snp.eps

hold on
chm = ras.chm10(jj,ii);
x = ras.x(ii);
y = ras.y(jj);
locs = locmax(chm,3,3);
[X,Y] = meshgrid(x,y);
hp = plot(X(locs),Y(locs),'.w','markersize',3);
print -depsc2 -r600 /Users/morsdorf/Graphiken/chm_snp_loc.eps