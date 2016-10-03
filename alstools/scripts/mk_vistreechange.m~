cc
cd /Users/morsdorf/Data/ForestGrowth/LFI_Wald
S = shaperead('LFI_Wald.shp');
cd /Users/morsdorf/Data/ForestGrowth

lfi = xlsread('/Users/morsdorf/data/ForestGrowth/ba_peri_20140127_v2.xlsx');
return
load trees2002
alltrees02 = trees;

load trees2010
alltrees10 = trees;
%[alltrees,statsall] = matchtrees(alltrees02,alltrees10,1,'als');
%save alltrees

load dsmdtm.mat
load alltrees
ras.dsm10(isnan(ras.dsm02)) = NaN;
ras.dtm10(isnan(ras.dsm02)) = NaN;

load polygon2
ypol1(1) = ypol1(1)-20;
ypol2(2) = ypol2(2)-20;

ii = inpolygon(ras.X(:),ras.Y(:),xpol1([1:3,1]),ypol1([1:3,1]));
ras.dtm10(ii) = NaN;
ii = inpolygon(ras.X(:),ras.Y(:),xpol2([1:3,1]),ypol2([1:3,1]));
ras.dtm10(ii) = NaN;
siz = 0.5;

%shapewrite(S2,'Polygon_SNP');

% clean up tree structure
ii = alltrees.h < 40 & alltrees.h > 3;
alltrees = subsetraw(alltrees,ii);
figure;
imagesc(ras.x,ras.y,ras.dsm10-ras.dtm10);hold on;caxis([0 30]);
colormap(gray);
plot(alltrees.x,alltrees.y,'.g','markersize',siz);
plot(alltrees.alsx,alltrees.alsy,'.r','markersize',siz);
plot(lfi(:,5),lfi(:,6),'xy','markersize',siz);
hl = line([alltrees.alsx,alltrees.x]',[alltrees.alsy,alltrees.y]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
axis([min(ras.x) max(ras.x) min(ras.y) max(ras.y)]);
axis equal;axis tight;axis xy;
mytitle('Matched tree locations - 2002 - 2010');
legend('2002','2010','location','northwest');
%print -depsc2 -r1800 treelocs_matchall.eps 

return
dh = alltrees.alsh - alltrees.h;
dh(abs(dh)>40) = NaN;
dh(dh>3) = 3;
dh(dh<-3) = -3;
figure;
ii = ~isnan(alltrees.h);
myscatter3(alltrees.x(ii),alltrees.y(ii),alltrees.h(ii),dh(ii),jet(256),3);
view(2);
swisstick
grid on
box on
colormap(flipud(jet(256)));
colorbar
mytitle(['Tree height change 2010 - 2002']);
axis equal
axis tight
orient landscape
wysiwyg
print -dtiff -r600 treeheightchangeall.tif

figure;
load lwf00
shadmod(ras.x,ras.y,ras.dtm10);
hold on
alltrees.z = interp2(ras.X,ras.Y,ras.dtm10,alltrees.x,alltrees.y);
ii = [S.X] >= nanmin(ras.X(:)) & [S.X] <= nanmax(ras.X(:)) & [S.Y] >= nanmin(ras.Y(:)) & [S.Y] <= nanmax(ras.Y(:));
SSNP.X = [S(ii).X];SSNP.Y = [S(ii).Y];
SSNP.Z = interp2(ras.X,ras.Y,ras.dtm10,SSNP.X,SSNP.Y);
hp = myscatter3(alltrees.x,alltrees.y,alltrees.z+3,dh,flipud(ocean2(256)),2);
colormap(flipud(ocean2(256)))
view(2);
swisstick
axis equal
axis tight
orient landscape
wysiwyg
hold on
plot3(SSNP.X,SSNP.Y,SSNP.Z,'xow');
dtm.x = ras.x;
dtm.y = ras.y;
dtm.z = ras.dtm10;
plotrect3d([min(lwf00.x) max(lwf00.x) min(lwf00.y) max(lwf00.y)],dtm);
hcb = colorbar;
axes(hcb);ylabel('ALS based tree height change [m]');
print -depsc2 -r600 treeheightchangeall3d_1.eps
colormap(gray);
print -depsc2 -r600 treeheightchangeall3d_2.eps


return
figure;
myimage(ras.x,ras.y,ras.dsm10-ras.dtm10);
colormap(jet);
caxis([0 40]);
colorbar;
hold on
plot(alltrees.x,alltrees.y,'.w','markersize',0.05);
print -depsc2 -r1200 treelocsgrowthall.eps
if 0 
    [slo,asp] = slaspect(ras.x,ras.y,ras.dtm10);
    save sloasp slo asp
else
    load sloasp;
    alltrees.slo = interp2(ras.X,ras.Y,slo,alltrees.x,alltrees.y);
    alltrees.asp = interp2(ras.X,ras.Y,asp,alltrees.x,alltrees.y);
end
load poly_ns;
figure;
subplot(3,1,1)
dh = alltrees.alsh - alltrees.h;
ii = abs(dh) < 40;
jj = inpolygon(alltrees.x,alltrees.y,xi,yi);
ii = ii & ~jj;
SLO = slo(ii);
DH = dh(ii);
ASP = asp(ii);
ALT = alltrees.z(ii);
alltrees.dh = dh;
trees = subsetraw(alltrees,ii);
res = 20;
sloi = 0-res/2:res:max(SLO)+res/2;
for i = 1:length(sloi)-1;
    ii = SLO >= sloi(i) & SLO < sloi(i+1);
    mdh(i) = nanmean(DH(ii));
    sdh(i) = nanstd(DH(ii));
end
%errorbar(0:res:max(SLO),mdh,sdh);
plot(0:res:max(SLO),mdh,'-x');
xlabel('Slope');
clear mdh sdh    

subplot(3,1,2)
res = 45;
aspi = 0-res/2:res:max(ASP)+res/2;
for i = 1:length(aspi)-1;
    ii = ASP >= aspi(i) & ASP < aspi(i+1);
    mdh(i) = nanmean(DH(ii));
    sdh(i) = nanstd(DH(ii));
end
%errorbar(0:res:max(ASP),mdh,sdh);
plot(0:res:max(ASP),mdh,'-x');
xlabel('Aspect');


subplot(3,1,3)
res = 50;
hghi = nanmin(alltrees.z)-res/2:res:nanmax(alltrees.z)+res/2;
for i = 1:length(hghi)-1;
    ii = ALT >= hghi(i) & ALT < hghi(i+1);
    mdh(i) = nanmean(DH(ii));
    sdh(i) = nanstd(DH(ii));
end
%errorbar(0:res:max(ASP),mdh,sdh);
plot(nanmin(alltrees.z):res:nanmax(alltrees.z),mdh,'-x');
xlabel('Altitude [m]');
print -depsc2 -r300 heightchangeslaspectalti.eps

% save trees to file

trees2shp('SingleTreesOfenpass2002_2010',trees);
figure;
%are = [813350 814400 171250 172000];
are = [813530 813670 171600 171780];
ii = alltrees.x >= are(1) & alltrees.x <= are(2) & alltrees.y >= are(3) & alltrees.y <= are(4);
trees = subsetraw(alltrees,ii);
ii = ras.x >= are(1) & ras.x <= are(2);
jj = ras.y >= are(3) & ras.y <= are(4);
clear dtm;dtm.x = ras.x(ii);dtm.y = ras.y(jj);dtm.z = ras.dtm10(jj,ii);
load lwf00;
dhs = trees.alsh - trees.h;
dhs(abs(dhs)>40) = NaN;
dhs(dhs>3) = 3;
dhs(dhs<-3) = -3;
visplot3(trees,dtm,dhs);
hold on 
plot3(lwf00.x,lwf00.y,lwf00.z+0.5,'.w','markersize',7);
print -depsc2 -r300 treeheightchangelwf3d.eps


figure;
dh = alltrees.alsh - alltrees.h;
ii = abs(dh) < 40;
H = alltrees.h(ii);DH = dh(ii);
res = 0.5;
dth = nanmin(H)-res/2:res:nanmax(H)+res/2;
for i = 1:length(dth)-1;
    ii = H >= dth(i) & H < dth(i+1);
    mdh(i) = nanmean(DH(ii));
    sdh(i) = nanstd(DH(ii));
end
errorbar(nanmin(H):res:nanmax(H),mdh,sdh);
mytitle('Tree height vs. \Delta height');
print -depsc2 -r600 treegrowthvsheight.eps