% script for producing tree height regression for NCCOOPS GCB paper

cc
% load WSL tree data

dat = csvread(['/Users/morsdorf/Teaching/MSc/MA_SW/WSL/Umgewandelte_Originale/', ...
              'Angepasste_Hoeheninventur_2005_woheader.csv']);
lwf.x = dat(:,13);
lwf.y = dat(:,14);
lwf.z = dat(:,15);
lwf.dbh = dat(:,3)/pi;
lwf.h = dat(:,5);

dat2 = csvread(['/Users/morsdorf/Teaching/MSc/MA_SW/WSL/Umgewandelte_Originale/', ...
              'Angepasste_Hoeheninventur_2010_woheader.csv']);
lwf2.x = dat2(:,13);
lwf2.y = dat2(:,14);
lwf2.z = dat2(:,15);
lwf2.dbh = dat2(:,3)/pi;
lwf2.h = dat2(:,5);

load /Users/morsdorf/Teaching/MSc/MA_SW/Arbeit/MA_Walker_2013_DVD_2/0_Originaldaten/RASTER2010/dsm2010.mat
load /Users/morsdorf/Teaching/MSc/MA_SW/Arbeit/MA_Walker_2013_DVD_2/0_Originaldaten/RASTER2010/dtm2010.mat

dx = 2;
ii = dsm.x >= min(lwf.x)-dx & dsm.x <= max(lwf.x)+dx;
jj = dsm.y >= min(lwf.y)-dx & dsm.y <= max(lwf.y)+dx;
dsm.x = dsm.x(ii);
dsm.y = dsm.y(jj);
dsm.z = dsm.z(jj,ii);

ii = dtm.x >= min(lwf.x)-dx & dtm.x <= max(lwf.x)+dx;
jj = dtm.y >= min(lwf.y)-dx & dtm.y <= max(lwf.y)+dx;
dtm.x = dtm.x(ii);
dtm.y = dtm.y(jj);
dtm.z = dtm.z(jj,ii);
dtm.z = inpaint_nans(dtm.z,4);
chm = dsm;
chm.z = dsm.z - dtm.z;
lwf = lwf2;
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);
% use only trees with height measurement

ii = lwf.h ~= -9999;
lwf = subsetraw(lwf,ii);
ii = lwf.dbh > 12;
lwf = subsetraw(lwf,ii);

lwf.als = interp2(chm.X,chm.Y,chm.z,lwf.x,lwf.y);
figure(1);clf
plot([0,18],[0,18],'-k','linewidth',2);
hold on; axis([0 18 0 18])
mytitle('Comparison of ALS canopy height at field measurement locations');
add_regress(lwf.h,lwf.als,{'LWF Field Measured Tree Height [m]','ALS Canopy Height [m]'});
print -depsc2 -r600 /Users/morsdorf/Graphiken/gcb_treeheight_img1.eps

figure(2);clf;myimage(chm);
ii = locmax(chm.z,3,0.7,3);
hold on
hp(1) = plot(chm.X(ii),chm.Y(ii),'.','color',[0.9 0.9 0.9])
hp(2) = plot(lwf.x,lwf.y,'.k')
legend(hp,'ALS tree locations','Field measured tree locations'); 
mytitle('Canopy height model and tree locations');
hcb = colorbar;
axes(hcb);ylabel('Canopy height [m]');

print -depsc2 -r600 /Users/morsdorf/Graphiken/gcb_treeheight_img2.eps
% matchtrees
als.x = chm.X(ii);
als.y = chm.Y(ii);
als.h = chm.z(ii);

[trees,stats,id] = matchtrees(lwf,als,1);

ii = ~isnan(trees.h(:,1));
figure(3);clf;
plot([0,18],[0,18],'-k','linewidth',2);
hold on; axis([0 18 0 18])
mytitle('Comparison of ALS tree height with field measured tree height');
add_regress(trees.h(ii,1),trees.h(ii,2),{'LWF Field Measured Tree Height [m]','ALS Measured Tree Height [m]'});
csvwrite('/Users/morsdorf/Desktop/tree_heights_SNP_subset.csv',[trees.h(ii,1),trees.h(ii,2)]);
print -depsc2 -r600 /Users/morsdorf/Graphiken/gcb_treeheight_img3.eps