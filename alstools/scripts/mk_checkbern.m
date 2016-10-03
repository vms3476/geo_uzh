cc
% script for visualizing testdata of berner jura
%cd /Volumes/Data2/BernerJura/DOM/Grid
%dom = xyz2ras('1125_14.xyz');
%cd /Volumes/Data2/BernerJura/DTM/Grid
%dtm = xyz2ras('1125_14.xyz');
%save /Volumes/Data2/BernerJura/DOMDTM.mat dom dtm
load /Volumes/Data2/BernerJura/DOMDTM.mat 
if 1
%make pointdensity maps
are = [2571875.3   2576249.8   1224000.3   1226999.8];
%[x,y,numechos,numshots] = pointdens_tiles('/Volumes/Data2/BernerJura/Punktwolke_klassiert',are,1);
%save /Volumes/Data2/BernerJura/Punktwolke_klassiert/pdens.mat x y numechos numshots
load /Volumes/Data2/BernerJura/Punktwolke_klassiert/pdens.mat

% figure 1 - shaded dtm with colormap of pointdens < 4 must be red
[X,Y] = meshgrid(x,y);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
newz = interp2(dtm.X,dtm.Y,dtm.z,X,Y);
figure(1);clf
subplot(2,1,1)
shadmod(x,y,newz,numshots);axis xy;axis equal;axis tight;swisstick;
title('Schussdichte Kachel Cortebert')
caxis([0 20])
hcb = colorbar;
set(hcb,'ylim',[0 18],'ytick',[0:2:18]);
axes(hcb);ylabel('Schussdichte [m^{-2}]');
subplot(2,1,2)
shadmod(x,y,newz,numechos);axis xy;axis equal;axis tight;swisstick;
title('Echodichte Kachel Cortebert')
caxis([0 20])
colormap([(hsv(10))])
hcb(2) = colorbar;
set(hcb,'ylim',[0 18],'ytick',[0:2:18]);
axes(hcb(2));ylabel('Echodichte [m^{-2}]');


print -r600 -depsc2 PunktdichteCortebert.eps

% figure 2 - histogramm der punktdichte
figure(2);clf
subplot(2,1,1)
[N,X] = hist(numshots(:),[0:0.5:30]);
N = N/length(numshots(:))*100;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 15]);
xlabel('Schuesse pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = numshots;
NH(isnan(numshots)) = [];
title(['Histogramm Schussdichte Kachel Cortebert, Mittel: ',num2str(mean(NH(:))), ...
       ', Median: ',num2str(median(NH(:)))]);
subplot(2,1,2)
[N,X] = hist(numechos(:),[0:0.5:30]);
N = N/length(numechos(:))*100;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 15]);
xlabel('Echos pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = numechos;
NH(isnan(numechos)) = [];
title(['Histogramm Echodichte Kachel Cortebert, Mittel: ',num2str(mean(NH(:))), ...
       ', Median: ',num2str(median(NH(:)))]);

% lage der referenzdaten
end

figure(3)
load /Volumes/Data2/BernerJura/Referenz/punkte.dat
buf = 5;
are = [min(punkte(:,2))-buf max(punkte(:,2))+buf min(punkte(:,3))-buf max(punkte(:,3))+buf];
%raw = getrawlas('/Volumes/Data2/BernerJura/Punktwolke_klassiert',are);
ii = find(dom.x >= are(1) & dom.x <= are(2));
jj = find(dom.y >= are(3) & dom.y <= are(4));

DOM.z = dom.z(jj,ii);
DOM.x = dom.x(ii);
DOM.y = dom.y(jj);
shadmod(DOM);
hold on
[DOM.X,DOM.Y] = meshgrid(DOM.x,DOM.y);
DTM.z = dtm.z(jj,ii);
DTM.x = dtm.x(ii);
DTM.y = dtm.y(jj);
[DTM.X,DTM.Y] = meshgrid(DTM.x,DTM.y);
pz = interp2(DOM.X,DOM.Y,DOM.z,punkte(:,2),punkte(:,3));
plot3(punkte(:,2),punkte(:,3),punkte(:,4),'ow');
plot3(punkte(:,2),punkte(:,3),punkte(:,4),'xk');
text(punkte(:,2),punkte(:,3),punkte(:,4),[num2str(punkte(:,1))]);
title('Haus in Kachel Cortebert');
swisstick
print -depsc2 -r600 DSM_Haus_Referenz_Cortebert_topview.eps
axis equal;axis tight;
view(14,40);
print -depsc2 -r600 DSM_Haus_Referenz_Cortebert_sideview.eps
figure(4);
load /Volumes/Data2/BernerJura/Punktwolke_klassiert/raw_haus.mat
mex = mean(raw.x);
mey = mean(raw.y)
myscatter3(raw.x-mex,raw.y-mey,raw.z,raw.z);
hold on
plot3(punkte(:,2)-mex,punkte(:,3)-mey,punkte(:,4),'ok');
axis equal;axis tight;
view(2)
print -depsc2 -r600 raw_Haus_Referenz_Cortebert_topview.eps
view(-5.5,0);
print -depsc2 -r600 raw_Haus_Referenz_Cortebert_sideview.eps
mat3d2osg('Haus_Cortebert')
diffp = SurveyPT_vs_pointcloud(punkte(:,2:4),[raw.x,raw.y,raw.z],1,2);
newp = punkte(diffp(:,2)<0.5,2:4);
raw.h = raw.z-interp2(DTM.X,DTM.Y,DTM.z,raw.x,raw.y);
ii = raw.h > 5;
[R, t, corr, error, data2] = icp2(newp',[raw.x(ii),raw.y(ii),raw.z(ii)]',0.5);
results = [punkte';diffp']';
fid = fopen('/Volumes/Data2/BernerJura/Punktwolke_klassiert/results_match.csv','w');
for i = 1:length(results)
    fprintf(fid,'%4.0f %8.2f %8.2f %4.2f %2.4f %2.4f\n',results(i,:));
end
fclose(fid);
xlswrite('/Volumes/Data2/BernerJura/Punktwolke_klassiert/results_match.xls',results);
