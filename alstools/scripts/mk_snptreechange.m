% script for producing tree height change for swiss national park
%load old modelss
cc
cd /Users/morsdorf/Data/ForestGrowth
if 0
    load ras.mat
    chm = dsm;
    chm.z = dsm.z - dtm.z;
    % load tree reference data
    dat = csvread(['/Users/morsdorf/Teaching/MSc/MA_SW/WSL/Umgewandelte_Originale/', ...
                   'Angepasste_Hoeheninventur_2005_woheader.csv']);
    dat(dat == -9999) = NaN;
    lwf05.x = dat(:,13);
    lwf05.y = dat(:,14);
    lwf05.z = dat(:,15);
    lwf05.dbh = dat(:,3)/pi;
    lwf05.h = dat(:,5);
    lwf05.id1 = dat(:,1);
    lwf05.id2 = dat(:,2);
    
    dat2 = csvread(['/Users/morsdorf/Teaching/MSc/MA_SW/WSL/Umgewandelte_Originale/', ...
                    'Angepasste_Hoeheninventur_2010_woheader.csv']);
    dat2(dat2 == -9999) = NaN;
    lwf10.x = dat2(:,13);
    lwf10.y = dat2(:,14);
    lwf10.z = dat2(:,15);
    lwf10.dbh = dat2(:,3)/pi;
    lwf10.h = dat2(:,5);
    lwf10.id1 = dat2(:,1);
    lwf10.id2 = dat2(:,2);
    
    dat3 = csvread(['/Users/morsdorf/Teaching/MSc/MA_SW/WSL/Umgewandelte_Originale/', ...
                    'Angepasste_Hoeheninventur_2000_woheader.csv']);
    dat3(dat3 == -9999) = NaN;
    lwf00.x = dat3(:,13);
    lwf00.y = dat3(:,14);
    lwf00.z = dat3(:,15);
    lwf00.dbh = dat3(:,3)/pi;
    lwf00.h = dat3(:,5);
    lwf00.id1 = dat3(:,1);
    lwf00.id2 = dat3(:,2);
    
    % find heights and dbh of corresponding trees
    lwf10.h00 = lwf10.x*NaN;
    lwf10.dbh00 = lwf10.x*NaN;
    lwf10.h05 = lwf10.x*NaN;
    lwf10.dbh05 = lwf10.x*NaN;
    for i = 1:length(lwf10.x);
        ii = find(lwf00.id1 == lwf10.id1(i));
        if ~isempty(ii)
            lwf10.h00(i) = lwf00.h(ii);
            lwf10.dbh00(i) = lwf00.dbh(ii);
        end
        ii = find(lwf05.id1 == lwf10.id1(i));
        if ~isempty(ii)
            lwf10.h05(i) = lwf05.h(ii);
            lwf10.dbh05(i) = lwf05.dbh(ii);
        end
    end
    ras.dtm10 = interp2(dtm.X,dtm.Y-0.5,dtm.z,ras.X,ras.Y,'cubic');
    ras.dsm10 = interp2(dsm.X,dsm.Y-0.5,dsm.z,ras.X,ras.Y,'cubic');
    %s = sysshift(ras.dtm02(750:1250,2000:3000),ras.dtm10(750:1250,2000:3000));
    ras.dsm10(isnan(ras.dsm02)) = NaN;
    ras.dtm10(isnan(ras.dsm02)) = NaN;
    ras.chm02 = ras.dsm02-ras.dtm02;
    ras.chm10 = ras.dsm10-ras.dtm10;
    ii = ras.x >= min(lwf10.x)-10 & ras.x <= max(lwf10.x)+10;
    jj = ras.y >= min(lwf10.y)-10 & ras.y <= max(lwf10.y)+10;
    chm02_1m.X = ras.X(jj,ii);
    chm02_1m.Y = ras.Y(jj,ii);
    chm02_1m.x = ras.x(ii);
    chm02_1m.y = ras.y(jj);
    chm02_1m.z = ras.chm02(jj,ii);
    
    ii = chm.x >= min(lwf10.x)-10 & chm.x <= max(lwf10.x)+10;
    jj = chm.y >= min(lwf10.y)-10 & chm.y <= max(lwf10.y)+10;
    
    chm10.X = chm.X(jj,ii);
    chm10.Y = chm.Y(jj,ii);
    chm10.x = chm.x(ii);
    chm10.y = chm.y(jj);
    chm10.z = chm.z(jj,ii);
    DSM = dsm;
    DTM = dtm;
    load LWF_DSM
    load LWF_DTM
    chm = dtm;
    chm.z = dsm.z(1:end-1,:) - dtm.z;
    ii = chm.x >= min(lwf10.x)-10 & chm.x <= max(lwf10.x)+10;
    jj = chm.y >= min(lwf10.y)-10 & chm.y <= max(lwf10.y)+10;
    [chm.X,chm.Y] = meshgrid(chm.x,chm.y);
    chm02.X = chm.X(jj,ii);
    chm02.Y = chm.Y(jj,ii);
    chm02.x = chm.x(ii);
    chm02.y = chm.y(jj);
    chm02.z = chm.z(jj,ii);
    clear ras;
    % load raw data

    raw10 = readlas('raw2010.las');
    RAW10.x = raw10.x;RAW10.y = raw10.y;RAW10.z = raw10.z;
    clear raw10;raw10 = RAW10;clear RAW10;
    load raw2002;raw02 = raw2002;clear raw2002;
    % load original raw data for subset

    are = [min(raw10.x) max(raw10.x) min(raw10.y) max(raw10.y)];
    raw02_old = raw02;
    clear raw02;raw02.x = [];
    raw02.x = [];
    raw02.y = [];
    raw02.z = [];
    raw02.f = [];
    fpath = '/Users/morsdorf/Data/ForestGrowth/rohdaten/first';
    filez = dir([fpath,'/*.ASC']);
    for i = 1:length(filez)
        dat = load([fpath,'/',filez(i).name]);
        ii = dat(:,2) >= are(1) & dat(:,2) <= are(2) & dat(:,1) >= are(3) & dat(:,1) <= are(4);
        raw02.x = [raw02.x;dat(ii,2)];
        raw02.y = [raw02.y;dat(ii,1)];
        raw02.z = [raw02.z;dat(ii,3)];
        raw02.f = [raw02.f;ones(size(dat(ii,1)))];
    end
    
    fpath = '/Users/morsdorf/Data/ForestGrowth/rohdaten/last';
    filez = dir([fpath,'/*.ASC']);
    for i = 1:length(filez)
        dat = load([fpath,'/',filez(i).name]);
        ii = dat(:,2) >= are(1) & dat(:,2) <= are(2) & dat(:,1) >= are(3) & dat(:,1) <= are(4);
        raw02.x = [raw02.x;dat(ii,2)];
        raw02.y = [raw02.y;dat(ii,1)];
        raw02.z = [raw02.z;dat(ii,3)];
        raw02.f = [raw02.f;zeros(size(dat(ii,1)))];
    end
    
    % interpolate terrain heights
    raw10.tz = interp2(DTM.X,DTM.Y,DTM.z,raw10.x,raw10.y);
    [dtm.X,dtm.Y]= meshgrid(dtm.x,dtm.y);
    raw02.tz = interp2(dtm.X,dtm.Y,dtm.z,raw02.x,raw02.y);
    save treegrowth.mat
end
load treegrowth
are =   [min(lwf10.x)-10 max(lwf10.x)+10 min(lwf10.y)-10 max(lwf10.y)+10];
% 1st figure, compute and display height changes due to field measurements and DBH
if 1
    figure(1);clf
    dh10 = lwf10.h-lwf10.h00;
    dh51 = lwf10.h05-lwf10.h00;
    dh52 = lwf10.h-lwf10.h05;
    
    ddbh10 = lwf10.dbh-lwf10.dbh00;
    ddbh51 = lwf10.dbh05-lwf10.dbh00;
    ddbh52 = lwf10.dbh-lwf10.dbh05;
    subplot(2,1,1)
    bins = -1:0.1:1;
    [fz,height] = hist(dh10,bins);
    [lz,height] = hist(dh51,bins);
    [mz,height] = hist(dh52,bins);
    cmap = gray(10);
    colormap(cmap([2,4,7],:));
    hb = bar(height,[fz;lz;mz]',1,'grouped');
    set(hb,'edgecolor','k')
    hl = legend(hb(1:3),['2010-2000, mean: ',num2str(round(nanmean(dh10)*100)/100), ...
                        ', std:',num2str(round(nanstd(dh10)*100)/100)], ...
                ['2005-2000, mean: ',num2str(round(nanmean(dh51)*100)/100), ...
                 ', std:',num2str(round(nanstd(dh51)*100)/100)], ...
                ['2010-2005, mean: ',num2str(round(nanmean(dh52)*100)/100), ...
                 ', std:',num2str(round(nanstd(dh52)*100)/100)], ...
                'location','northwest');
    xlabel('Tree height change [m]');
    mytitle('Tree height change distribution for three periods');
    ylabel('Number of trees');
    set(gca,'xlim',[min(bins) max(bins)]);
    hlt = findobj(hl,'type','text');
    set(hlt,'fontsize',10);
    grid on
    
    subplot(2,1,2)
    bins = -1:0.1:1;
    [fz,height] = hist(ddbh10,bins);
    [lz,height] = hist(ddbh51,bins);
    [mz,height] = hist(ddbh52,bins);
    colormap(cmap([2,4,7],:));
    hb = bar(height,[fz;lz;mz]',1,'grouped');
    set(hb,'edgecolor','k')
    hl = legend(hb(1:3),['2010-2000, mean: ',num2str(round(nanmean(ddbh10)*100)/100), ...
                        ', std:',num2str(round(nanstd(ddbh10)*100)/100)], ...
                ['2005-2000, mean: ',num2str(round(nanmean(ddbh51)*100)/100), ...
                 ', std:',num2str(round(nanstd(ddbh51)*100)/100)], ...
                ['2010-2005, mean: ',num2str(round(nanmean(ddbh52)*100)/100), ...
                 ', std:',num2str(round(nanstd(ddbh52)*100)/100)], ...
                'location','northwest');
    xlabel('Tree diameter change [cm]');
    mytitle('Tree diameter change distribution for three periods');
    ylabel('Number of trees');
    set(gca,'xlim',[min(bins) max(bins)]);
    hlt = findobj(hl,'type','text');
    set(hlt,'fontsize',10);
    grid on
    print -depsc2 lwf_treegrowth.eps
end




% shift lwf trees in 2002 with 0.5 meters to the east
lwf00.x = lwf00.x + 0.5;
tree02 = locmax(chm02.z,3,0.45,3);
%tree02 = locmax(chm02.z,3,0.15,3);
disp([num2str(sum(tree02(:))),' Trees found in 2002'])
tree10 = locmax(chm10.z,5,0.75,3);
%tree10 = locmax(chm10.z,3,0.5,3);
disp([num2str(sum(tree10(:))),' Trees found in 2010'])

trees02.x = chm02.X(tree02);trees02.y = chm02.Y(tree02);trees02.h = chm02.z(tree02);
trees10.x = chm10.X(tree10);trees10.y = chm10.Y(tree10);trees10.h = chm10.z(tree10);
raw10.oz = raw10.z;
raw10.z = raw10.z - raw10.tz;
ii = raw10.z < 50;
raw10 = subsetraw(raw10,ii);
raw02.oz = raw02.z;
raw02.z = raw02.z - raw02.tz;
ii = raw02.z < 50;
raw02 = subsetraw(raw02,ii);
% compute avg. point densities
pd02 = length(raw02.x)/((max(raw02.x)-min(raw02.x))*(max(raw02.y)-min(raw02.y)));
pd10 = length(raw10.x)/((max(raw10.x)-min(raw10.x))*(max(raw10.y)-min(raw10.y)));
% make cluster analysis
if 0
    [X10,Y10,Z10] = segtree(raw10,[trees10.x,trees10.y,trees10.h],3,4,'gmm');
    TREES10 = crdata(X10,Y10,Z10);
    [X02,Y02,Z02] = segtree(raw02,[trees02.x,trees02.y,trees02.h],3,4,'gmm');
    TREES02 = crdata(X02,Y02,Z02);
    save TREES_GMM TREES02 TREES10 X10 Y10 Z10 X02 Y02 Z02
else
    load TREES_GMM
end
% compute point density for each tree segment
TREES10.pd = TREES10.x*NaN;
TREES02.pd = TREES02.x*NaN;
for i = 1:length(TREES10.x)
    if sum(isnan(TREES10.olx{i})) == 0 
        [K,A] = convhull(TREES10.olx{i},TREES10.oly{i});
        TREES10.pd(i) = TREES10.points(i)/A;
    end
end

for i = 1:length(TREES02.x)
    if sum(isnan(TREES02.olx{i})) == 0 
        [K,A] = convhull(TREES02.olx{i},TREES02.oly{i});
        TREES02.pd(i) = TREES02.points(i)/A;
    end
end

%TREES02.h = TREES02.h90;
%TREES10.h = TREES10.h90;
figure(2);clf
orient landscape
wysiwyg
ii = lwf00.dbh > 12;
jj = lwf10.dbh > 12;
siz = 8;
subplot(1,2,1)
imagesc(chm02.x,chm02.y,chm02.z);hold on
plot(lwf00.x(ii),lwf00.y(ii),'.g','markersize',siz);
plot(trees02.x,trees02.y,'.r','markersize',siz);
swisstick
axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis tight;axis xy;
mytitle('2002');
legend('Field','ALS','location','northwest');
subplot(1,2,2)
imagesc(chm10.x,chm10.y,chm10.z);hold on
plot(lwf10.x(ii),lwf10.y(ii),'.g','markersize',siz);
plot(trees10.x,trees10.y,'.r','markersize',siz);

swisstick
axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis tight;axis xy;
colormap(gray)
mytitle('2010');
print -depsc2 -r600 treeloc2.eps



LWF00 = lwf00;
LWF10 = lwf10;
% use only trees with DBH > 12 cm !!!!
lwf00 = subsetraw(lwf00,ii);
lwf10 = subsetraw(lwf10,jj);
% use only dominant trees (within 1 m);
lwf00 = idtreeclust(lwf00,1);
lwf10 = idtreeclust(lwf10,1);
ii = lwf00.ID == 1;
lwf00 = subsetraw(lwf00,ii);
ii = lwf10.ID == 1;
lwf10 = subsetraw(lwf10,ii);
lwf10.h(isnan(lwf10.h)) = 0;
lwf00.h(isnan(lwf00.h)) = 0;



[lwf00c,stats00c] = matchtrees(lwf00,trees02,1);
[lwf10c,stats10c] = matchtrees(lwf10,trees10,1);
[lwf00p,stats00p] = matchtrees(lwf00,TREES02,1);
[lwf10p,stats10p] = matchtrees(lwf10,TREES10,1);

if 0
    res = [1:2:15];
    sig = [0.1:0.05:2];
    stats00 = ones(length(res),length(sig),4);
    stats10 = ones(length(res),length(sig),4);
    
    for i = 1:length(res)
        for j = 1:length(sig)
            tree02 = locmax(chm02.z,res(i),sig(j),3);
            disp([num2str(sum(tree02(:))),' Trees found in 2002'])
            tree10 = locmax(chm10.z,res(i),sig(j),3);
            disp([num2str(sum(tree10(:))),' Trees found in 2010'])
            trees02.x = chm02.X(tree02);trees02.y = chm02.Y(tree02);trees02.h = chm02.z(tree02);
            trees10.x = chm10.X(tree10);trees10.y = chm10.Y(tree10);trees10.h = chm10.z(tree10);
            [dum,stats00c] = matchtrees(lwf00,trees02,1);
            stats00(i,j,1) = stats00c.dr;
            stats00(i,j,2) = stats00c.fp;
            stats00(i,j,3) = stats00c.fn;
            stats00(i,j,4) = stats00c.dh;
            [dum,stats10c] = matchtrees(lwf10,trees10,1);
            stats10(i,j,1) = stats10c.dr;
            stats10(i,j,2) = stats10c.fp;
            stats10(i,j,3) = stats10c.fn;
            stats10(i,j,4) = stats10c.dh;
        end
    end
    figure(11);clf
    for i = 1:4;subplot(2,2,i);imagesc(sig,res,squeeze(stats00(:,:,i)));
        ylabel('Kernel resolution');xlabel('Kernel shape');
        mytitle(str{i});colorbar;
    end
    print -depsc2 paraset_chm00.eps
    
    figure(12);clf
    for i = 1:4;subplot(2,2,i);imagesc(sig,res,squeeze(stats00(:,:,i)));
        ylabel('Kernel resolution');xlabel('Kernel shape');
        mytitle(str{i});colorbar;
    end
    print -depsc2 paraset_chm10.eps
end



figure(3);clf
orient landscape
wysiwyg
ii = ~isnan(lwf00c.alsx);
subplot(1,2,1)
imagesc(chm02.x,chm02.y,chm02.z);hold on
plot(lwf00c.x,lwf00c.y,'.g','markersize',siz);
plot(lwf00c.alsx(ii),lwf00c.alsy(ii),'.r','markersize',siz);
hl = line([lwf00c.alsx(ii),lwf00c.x(ii)]',[lwf00c.alsy(ii),lwf00c.y(ii)]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis tight;axis xy;
mytitle('2002');
legend('Field','ALS','location','northwest');
subplot(1,2,2)
ii = ~isnan(lwf10c.alsx);
imagesc(chm10.x,chm10.y,chm10.z);hold on
plot(lwf10c.x,lwf10c.y,'.g','markersize',siz);
plot(lwf10c.alsx(ii),lwf10c.alsy(ii),'.r','markersize',siz);
hl = line([lwf10c.alsx(ii),lwf10c.x(ii)]',[lwf10c.alsy(ii),lwf10c.y(ii)]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis tight;axis xy;
colormap(gray)
mytitle('2010');
hcb = colorbar('East');
axes(hcb);ylabel('Canopy Height [m]');set(gca,'yaxislocation','left','color','w','ycolor','w','xcolor','w');
set(gcf,'inverthardcopy','off','color','w');
setstrings(12)
print -depsc2 -r600 treeloc_match_chm.eps

figure(4);clf
orient landscape
wysiwyg
ii = ~isnan(lwf00p.alsx);
subplot(1,2,1)
imagesc(chm02.x,chm02.y,chm02.z);hold on
%plot(trees02.x,trees02.y,'.y','markersize',siz);
plot(lwf00p.x,lwf00p.y,'.g','markersize',siz);
plot(lwf00p.alsx(ii),lwf00p.alsy(ii),'.r','markersize',siz);
hl = line([lwf00p.alsx(ii),lwf00p.x(ii)]',[lwf00p.alsy(ii),lwf00p.y(ii)]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
%axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis xy;axis tight
%axis([813540 813630 171605 171660])
mytitle('2002');
legend('Field','ALS','location','northwest');

subplot(1,2,2)
ii = ~isnan(lwf10p.alsx);
imagesc(chm10.x,chm10.y,chm10.z);hold on
plot(lwf10p.x,lwf10p.y,'.g','markersize',siz);
plot(lwf10p.alsx(ii),lwf10p.alsy(ii),'.r','markersize',siz);
hl = line([lwf10p.alsx(ii),lwf10p.x(ii)]',[lwf10p.alsy(ii),lwf10p.y(ii)]');
set(hl,'linewidth',siz/10,'color','w');
swisstick
%axis([min(chm02.x) max(chm02.x) min(chm02.y) max(chm02.y)]);
axis equal;axis xy;axis tight;
%axis([813540 813630 171605 171660])
colormap(gray)
mytitle('2010');
hcb = colorbar('East');
axes(hcb);ylabel('Canopy Height [m]');set(gca,'yaxislocation','left','color','w','xcolor','w','ycolor','w');
setstrings(12)
set(gcf,'inverthardcopy','off','color','w');
print -depsc2 -r600 treeloc_match_plc.eps



% make regression of tree heights for each year
editfit(3,'mylin','a*x+b; a=1; b=1;'); 

figure(5);clf
subplot(1,2,1)
ii = ~isnan(lwf00c.alsh) & lwf00c.h < 30 & abs(lwf00c.alsh-lwf00c.h) < 5;

mytitle('2002');hold on
x = lwf00c.h(ii);
y = lwf00c.alsh(ii);
%[sf,gd] = fit( x, y, 'poly1');
%plot(sf, x, y)
ax = [0 20 0 20];
%plot(lwf00c.h(ii),lwf00c.alsh(ii),'.k');
%ezfit('mylin')
%showfit('mylin');
str = {'Field measured tree height [m]','ALS measured tree height, CHM [m]'};
[modchm02,gof02chm] = myregress(x,y,str,[1 NaN],ax);
%xlabel(str{1})
%ylabel(str{2})
subplot(1,2,2)
ii = ~isnan(lwf10c.alsh) & lwf10c.h < 30 & abs(lwf10c.alsh-lwf10c.h) < 5;
mytitle('2010');hold on
x = lwf10c.h(ii);
y = lwf10c.alsh(ii);
%[sf,gd] = fit( x, y, 'poly1');
%plot(sf, x, y)
%plot(lwf10c.h(ii),lwf10c.alsh(ii),'.k');
%ezfit('mylin')
%showfit('mylin');
str = {'Field measured tree height [m]','ALS measured tree height, CHM [m]'};
[modchm10,gof10chm] = myregress(x,y,str,[1 NaN],ax);
%xlabel(str{1})
%ylabel(str{2})
print -depsc2 -r600 treereg2_chm.eps




figure(6);clf
subplot(1,2,1)
ii = ~isnan(lwf00p.alsh) & abs(lwf00p.alsh-lwf00p.h) < 5;
mytitle('2002');hold on
[modplc02,gof02plc] = myregress(lwf00p.h(ii),lwf00p.alsh(ii),{'Field measured tree height [m]','ALS measured tree height, point cloud [m]'},[1 NaN],ax);

subplot(1,2,2)
ii = ~isnan(lwf10p.alsh) & abs(lwf10p.alsh-lwf10p.h) < 5;
mytitle('2010');hold on
[modplc10,gof10plc] = myregress(lwf10p.h(ii),lwf10p.alsh(ii),{'Field measured tree height [m]','ALS measured tree height, point cloud [m]'},[1 NaN],ax);

print -depsc2 -r600 treereg2_plc.eps
clear fldh alsh



% compute growth signal for field trees
k = 0;
for i = 1:length(lwf00.x);
    ii = find(lwf00.id1(i) == lwf10.id1);
    if ~isempty(ii)
        k = k + 1; 
        dlwf.x(k) = lwf10.x(ii);
        dlwf.y(k) = lwf10.y(ii);
        dlwf.h(k) = lwf10.h(ii);
        dlwf.dh(k) = lwf10.h(ii)-lwf00.h(i);
        dlwf.ddbh(k) = lwf10.dbh(ii)-lwf00.dbh(i);
        dlwf.dbh(k) = lwf10.dbh(ii);
    end
end
%dchm = matchtrees(trees02,trees10);
%dplc = matchtrees(TREES02,TREES10);
%dplc.h = dplc.alsh - dplc.h;
%dchm.h = dchm.alsh - dchm.h;

%chng_plc = matchtrees(dlwf,dplc);
%chng_chm = matchtrees(dlwf,dchm);



k = 0;
% link field trees for 2002 and 2010 (and thus ALS trees as well) % point cloud 
for i = 1:length(lwf10p.x);
    if ~isnan(lwf10p.alsh(i))
        k = k + 1;
        ii = find(lwf00p.id1 == lwf10p.id1(i));
        fldh(k,:) = [lwf00p.h(ii),lwf10p.h(i)];
        flddbh(k,:) = [lwf00p.dbh(ii),lwf10p.dbh(i)];
        alsh(k,:) = [lwf00p.alsh(ii),lwf10p.alsh(i)];
    end
end
% link field trees for 2002 and 2010 (and thus ALS trees as well) % chm
clear cfldh calsh
k = 0;
for i = 1:length(lwf10c.x);
    if ~isnan(lwf10c.alsh(i))
        k = k + 1;
        ii = find(lwf00c.id1 == lwf10c.id1(i));
        cfldh(k,:) = [lwf00c.h(ii),lwf10c.h(i)];
        cflddbh(k,:) = [lwf00c.dbh(ii),lwf10c.dbh(i)];
        calsh(k,:) = [lwf00c.alsh(ii),lwf10c.alsh(i)];
    end
end

figure(7);clf
orient landscape
wysiwyg
ax = [-3 3 -3 3];
subplot(1,2,2)
dhfld = fldh(:,2)-fldh(:,1);
dhals = alsh(:,2)-alsh(:,1);
%dhfld = chng_plc.dh(:);
%dhals = chng_plc.alsh(:);
ii = ~isnan(dhfld) & ~isnan(dhals) & abs(dhals) < 3 & dhfld >= -1.5;
myregress(dhfld(ii),dhals(ii),{'Field tree height change [m]','ALS tree height change, point cloud [m]'},[1 NaN], ax);
mytitle('\Delta h point cloud vs. \Delta h field');
%print -depsc2 -r600 treegrowth_reg_plc.eps

%figure(8);clf
subplot(1,2,1)
cdhfld = cfldh(:,2)-cfldh(:,1);
cdhals = calsh(:,2)-calsh(:,1);
%dhfld = chng_chm.dh(:);
%dhals = chng_chm.alsh(:);

nanmean(dhals);
ii = ~isnan(cdhfld) & ~isnan(cdhals) & abs(cdhals) < 3 & cdhfld >= -1.5;
myregress(cdhfld(ii),cdhals(ii),{'Field tree height change [m]','ALS tree height change, CHM [m]'},[1 NaN],ax);
%print -depsc2 -r600 treegrowth_reg_chm.eps
mytitle('\Delta h CHM vs. \Delta h field');
setstrings(12)
print -depsc2 -r600 treegrowth_2regs.eps
figure(9);clf
ddbhfld = flddbh(:,2)-flddbh(:,1);
dhals = alsh(:,2)-alsh(:,1);

%ddbhfld = chng_plc.ddbh(:);
%dhals = chng_plc.alsh(:);
nanmean(dhals);
ax = [-2 2 -4 4];
subplot(2,1,1)
ii = ~isnan(ddbhfld) & ~isnan(dhfld) & abs(dhfld) < 3;
myregress(ddbhfld(ii),dhfld(ii),{'Field DBH change [cm]','Field tree height change [m]'},[NaN,NaN],ax);
mytitle('\Delta h field vs. \Delta DBH field');
setstrings(12)

subplot(2,1,2)
ii = ~isnan(ddbhfld) & ~isnan(dhals) & abs(dhals) < 3;
myregress(ddbhfld(ii),dhals(ii),{'Field DBH change [cm]','ALS tree height change, point cloud [m]'},[NaN,NaN],ax);
mytitle('\Delta h point cloud vs. \Delta DBH field');
setstrings(12)

print -depsc2 -r600 treegrowth_regdbh.eps


figure
ax = [12 32 0 27];
subplot(2,2,3)
ii = ~isnan(lwf00p.alsh);
myregress(lwf00p.dbh(ii),lwf00p.alsh(ii),{'Field DBH','ALS derived tree height'},[NaN,NaN],ax,'power1');
mytitle('2000/2002');

subplot(2,2,4)
ii = ~isnan(lwf10p.alsh);
myregress(lwf10p.dbh(ii),lwf10p.alsh(ii),{'Field DBH','ALS derived tree height'},[NaN,NaN],ax,'power1');
mytitle('2010');


subplot(2,2,1)
ii = lwf00p.h ~= 0;
myregress(lwf00p.dbh(ii),lwf00p.h(ii),{'Field DBH','Field tree height'},[NaN,NaN],ax,'power1');
mytitle('2000/2002');
subplot(2,2,2)
ii = lwf10p.h ~= 0;
myregress(lwf10p.dbh(ii),lwf10p.h(ii),{'Field DBH','Field tree height'},[NaN,NaN],ax,'power1');
mytitle('2010');

print -depsc2 -r600 four_allometries.eps

figure;
str = {'2002','2010'}
for j = 1:2
    if j == 1;
        dh = lwf10p.h - lwf10p.alsh;
        pd = lwf10p.alspd;
    elseif j == 2
        dh = lwf00p.h - lwf00p.alsh;
        pd = lwf00p.alspd;
    end

    subplot(2,1,j)
    res = 1;
    dpd = 0:res:50;
    mh = ones(size(dpd))*NaN;
    sh = mh;
    for i = 1:length(dpd)
        ii = pd >= dpd(i)-res/2 & pd < dpd(i)+res/2;
        mh(i) = nanmean(dh(ii));
        sd(i) = nanstd(dh(ii));
    end
    errorbar(dpd,mh,sd); grid on;
    xlabel('Point density [m^{-2}]');
    ylabel('Tree height difference between field ans ALS trees');
    title(str{j});
end
print -depsc2 -r600 dh_vs_pointdensity.eps


% plot als based tree locs
figure(13)
plot(trees02.x,trees02.y,'xr',trees10.x,trees10.y,'og');
grid on
axis equal
axis tight
swisstick
legend('2002','2010');
mytitle('ALS based tree locations - CHM');
print -depsc2 als_treelocs_chm.eps

figure(14)
plot(TREES02.x,TREES02.y,'xr',TREES10.x,TREES10.y,'og');
grid on
axis equal
axis tight
swisstick
legend('2002','2010');
mytitle('ALS based tree locations - Point Cloud');
print -depsc2 als_treelocs_plc.eps

figure(15)
plot(trees02.x,trees02.y,'xr',TREES02.x,TREES02.y,'og');
grid on
axis equal
axis tight
swisstick
legend('CHM','PLC');
mytitle('ALS based tree locations - CHM/PLC 2002');
print -depsc2 als_treelocs_chmplc2002.eps

figure(16)
plot(trees10.x,trees10.y,'xr',TREES10.x,TREES10.y,'og');
grid on
axis equal
axis tight
swisstick
legend('CHM','PLC');
mytitle('ALS based tree locations - CHM/PLC 2010');
print -depsc2 als_treelocs_chmplc2010.eps

% match ALS based tree heights
dtrees = matchtrees(trees02,trees10,1);
dh = dtrees.alsh - dtrees.h;
figure(9);clf
ii = ~isnan(dh);
hist(dh(ii),100);

figure(10)
subplot(2,1,1);
diffh = lwf00p.alsh-lwf00c.alsh;
hist(diffh,-2:0.2:5)
set(gca,'xlim',[-2,5]);
xlabel('Height difference point cloud to CHM [m]');
ylabel('Number of trees')
mytitle(['2002 - mean difference: ',num2str(nanmean(diffh))])

subplot(2,1,2);
diffh = lwf10p.alsh-lwf10c.alsh;
hist(diffh,-2:0.2:5)
xlabel('Height difference point cloud to CHM [m]');
ylabel('Number of trees')
set(gca,'xlim',[-2,5]);
mytitle(['2010 - mean difference: ',num2str(nanmean(diffh))])
print -depsc2 -r600 treeheight_CHMvsPLC.eps
