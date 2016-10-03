cc
are = [669100 669100+1400 258800 259300];
are = [668900 675500 258300 260200];
load /Volumes/Data1/Laegeren/Spring_2010/DTM/dtm.mat
ii = dtm.x >= are(1) & dtm.x <= are(2);
jj = dtm.y >= are(3) & dtm.y <= are(4);
dtm.x = dtm.x(ii);
dtm.y = dtm.y(jj);
dtm.z = dtm.z(jj,ii);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
load /Volumes/Data1/Laegeren/Summer_2010/dsm_fe/dsm.mat
ii = dsm.x >= are(1) & dsm.x <= are(2);
jj = dsm.y >= are(3) & dsm.y <= are(4);
dsm.x = dsm.x(ii);
dsm.y = dsm.y(jj);
dsm.z = dsm.z(jj,ii);
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
dtmz = interp2(dtm.X,dtm.Y,dtm.z,dsm.X,dsm.Y);
chm = dsm;
chm.z = dsm.z - dtmz;
if 1
figure(1)
imagesc(chm.x,chm.y,chm.z);
axis equal;axis tight;axis xy;
swisstick;
colormap([1,1,1;myspecmap]);
caxis([0,50])
hcb = colorbar;
axes(hcb);
ylabel('Height [m]');
set(hcb,'ytick',[0:5:50]);
set(gcf,'renderer','zbuffer');
orient landscape
wysiwyg
print -depsc2 -r900 /Users/morsdorf/Desktop/laegeren_chm.eps
end
figure(3)
shadmod(dtm)
swisstick;box on
colormap(terrain);
hcb = colorbar;
axes(hcb);
ylabel('Height [m]');
set(gcf,'renderer','zbuffer');
orient landscape
wysiwyg
print -depsc2 -r900 /Users/morsdorf/Desktop/laegeren_dtm.eps
dsm.z = interp2(dsm.X,dsm.Y,dsm.z,dtm.X,dtm.Y);
figure(2)
ii = isnan(dtm.z);
dsm.z(ii) = NaN;
dtm.z = double(dsm.z);
shadmod(dtm)
swisstick
box on
colormap(terrain);
hcb = colorbar;
axes(hcb);
ylabel('Height [m]');
set(gcf,'renderer','zbuffer');
orient landscape
wysiwyg
print -depsc2 -r900 /Users/morsdorf/Desktop/laegeren_dsm.eps

