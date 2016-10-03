%make canopy height model of SNP as of 2010

load /Users/morsdorf/Desktop/ForestGrowth/dsm2010
load /Users/morsdorf/Desktop/ForestGrowth/dtm2010
are = [812800 815000 171000 172500];

ii = dtm.x >= are(1) & dtm.x <= are(2);
jj = dtm.y >= are(3) & dtm.y <= are(4);
dtm.x = dtm.x(ii);
dtm.y = dtm.y(jj);
dtm.z = dtm.z(jj,ii);


ii = dsm.x >= are(1) & dsm.x <= are(2);
jj = dsm.y >= are(3) & dsm.y <= are(4);
dsm.x = dsm.x(ii);
dsm.y = dsm.y(jj);
dsm.z = dsm.z(jj,ii);
dtm.z = inpaint_nans(dtm.z,4);
dsm.z = inpaint_nans(dsm.z,4);
chm.z = dsm.z-dtm.z;
myimage(dsm.x,dsm.y,chm.z);
caxis([0 30]);colormap([1,1,1;ocean2(15)])
title('Canopy Height Model SNP - Stabelchod');
hcb = colorbar;
yli = get(hcb,'ylabel')
set(yli,'string','Height [m]');
wysiwyg;
swisstick;
grid on;
print -depsc2 -r1200 /Users/morsdorf/Desktop/chm_snp_2010.eps
return
[X,Y] = meshgrid(dsm.x,dsm.y);

ii = locmax(chm.z,1,3);
tx = X(ii)+0.5;
ty = Y(ii)+0.5;
th = chm.z(ii);
hold on
plot(tx,ty,'.k','markersize',1);
print -depsc2 -r1200 /Users/morsdorf/Desktop/chm_snp_2010_seed.eps

%load /Volumes/Data1/KlostersDavosSNP/PUNKTWOLKE/lwf.mat
%tz = interp2(X,Y,dtm.z,raw.x,raw.y);
%raw.z = raw.z - tz;
%raw. z = double(raw.z);
%chm.x = dsm.x;
%chm.y = dsm.y;
%chm.X = X;
%chm.Y = Y;
%trees = segtree2(chm,raw,1,3,4);