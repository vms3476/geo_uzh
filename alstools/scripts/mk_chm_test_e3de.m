cc
are = [669100 669100+1400 258800 259300];
are = [668900 675500 258300 260200];
are = [669600      669900      258900      259200];
load /Volumes/Data1/Laegeren/Spring_2010/DTM/dtm.mat
trees = load('Users/morsdorf/Desktop/trees.csv');
ii = trees(:,1) >= are(1) & trees(:,1) <= are(2) & trees(:,2) >= are(3) & trees(:,2) ...
     <= are(4);
rtrees = trees(ii,:);
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
dtm.z = interp2(dtm.X,dtm.Y,dtm.z,dsm.X,dsm.Y);
chm = dsm;
chm.z = dsm.z - dtm.z;
imagesc(chm.x,chm.y,chm.z);
axis equal;axis tight;axis xy;
swisstick('german');
colormap([1,1,1;ocean2]);
caxis([0,50])
hold on
for i = 1:length(rtrees)
    circle(rtrees(i,4),rtrees(i,1),rtrees(i,2),'r',1.1,20);
end
load /Users/morsdorf/Desktop/trees_laegeren_own.mat
for i = 1:length(trees.x)
    circle(trees.dia(i),trees.x(i),trees.y(i),'b',1.1,20);
end
title('Vegetationshoehenmodell Laegeren');
hcb = colorbar;
axes(hcb);
ylabel('Hoehe [m]');
set(hcb,'ytick',[0:5:50]);



return

%set(gcf,'renderer','painters');
%orient portrait
%print -depsc2 -r600 /Users/morsdorf/Desktop/young_laegeren.eps