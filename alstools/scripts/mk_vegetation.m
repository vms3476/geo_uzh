% script for producing mosaics of fcover and LAI for Bernese Jura 650 km2!!!
close all;clear all
% make maps of products
res = 3;
if 0
    cd /Volumes/Data3/BernerJura/Punktwolke_klassiert/
    filez = dir('/Volumes/Data3/BernerJura/Punktwolke_klassiert/*_veg.mat');
    x = 2555000:res:2621000;
    y = 1211000:res:1245000;
    [X,Y] = meshgrid(x,y);
    FCO = ones(size(X))*NaN;
    LAI = FCO;
    TER = FCO;
    VEH = FCO;
    PDM = FCO;
    FCT = FCO;
    GPS = FCO;
    for i = 1:length(filez)
        disp(filez(i).name);
        load(filez(i).name);
        II = find(x >= min(xg)-res*0  &  x <= max(xg)+res*0);
        JJ = find(y >= min(yg)-res*0  &  y <= max(yg)+res*0);
        [XG,YG] = meshgrid(xg,yg);
        [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
        [m,n,o] = size(veg);
        if res <= 5
            if min(size(XG)) > 1 & o == 7 
                PDM(JJ,II) = interp2(XG,YG,veg(:,:,1),X(JJ,II),Y(JJ,II),'*cubic');
                FCT(JJ,II) = interp2(XG,YG,veg(:,:,2),X(JJ,II),Y(JJ,II),'*cubic');
                FCO(JJ,II) = interp2(XG,YG,veg(:,:,3),X(JJ,II),Y(JJ,II),'*cubic');
                LAI(JJ,II) = interp2(XG,YG,veg(:,:,4),X(JJ,II),Y(JJ,II),'*cubic');                
                VEH(JJ,II) = interp2(XG,YG,veg(:,:,6),X(JJ,II),Y(JJ,II),'*cubic');
                GPS(JJ,II) = interp2(XG,YG,veg(:,:,7),X(JJ,II),Y(JJ,II),'*cubic');
                II = find(x >= min(dtm.x)-res*0  &  x <= max(dtm.x)+res*0);
                JJ = find(y >= min(dtm.y)-res*0  &  y <= max(dtm.y)+res*0);
                TER(JJ,II) = interp2(dtm.X,dtm.Y,dtm.z,X(JJ,II),Y(JJ,II),'*cubic',-9999);
            end
        else
            for k = 1:length(II)
                for l = 1:length(JJ)
                    ii = find(xg >= x(II(k))-(res/2)  &  xg <= x(II(k))+(res/2));
                    jj = find(yg >= y(JJ(l))-(res/2)  &  yg <= y(JJ(l))+(res/2));
                    PDM(JJ(l),II(k)) = nanmean(nanmean(veg(jj,ii,1)));
                    FCO(JJ(l),II(k)) = nanmean(nanmean(veg(jj,ii,3)));
                    LAI(JJ(l),II(k)) = nanmean(nanmean(veg(jj,ii,4)));
                    VEH(JJ(l),II(k)) = nanmean(nanmean(veg(jj,ii,6)));
                    TER(JJ(l),II(k)) = nanmean(nanmean(dtm.z(jj,ii)));
                end
            end
        end
    end
    eval(['save /Volumes/Data3/BernerJura/vegetation_res',num2str(res),'.mat x y FCT FCO LAI VEH TER GPS PDM'])
end
if 1
    load ~/Projects/Archived/BernJura/forestclasses_3m.mat 
    eval(['load ~/Projects/Archived/BernJura/vegetation_res',num2str(res),'.mat'])
    dat = load('~/Projects/Archived/BernJura/flighttimes.txt');
    load /Users/morsdorf/Desktop/MODIS_LAI/TerraAcqua4Days/modis_lai2.mat
    mlai = lai; clear lai;
    [mlai.X,mlai.Y] = meshgrid(mlai.x{1},mlai.y{1});
    [cls.X,cls.Y] = meshgrid(cls.x+2000000,cls.y+1000000);
    [X,Y] = meshgrid(x,y);
    CLS = interp2(cls.X,cls.Y,cls.dat,X,Y);
    mlai.dat{1}(mlai.dat{1} >= 249 & qc{1} > 0) = NaN;
    MLAI = interp2(mlai.X,mlai.Y,[double(mlai.dat{1})],X,Y);
    mlai.dat{2}(mlai.dat{2} >= 249 & qc{2} > 0) = NaN;
    MLAI2 = interp2(mlai.X,mlai.Y,[double(mlai.dat{2})],X,Y);
    ii = mlai.x{1} >= min(x) & mlai.x{1} <= max(x);
    jj = mlai.y{1} >= min(y) & mlai.y{1} <= max(y);
    mlai.X = mlai.X(jj,ii);
    mlai.Y = mlai.Y(jj,ii);
    for i = 1:length(mlai.x);
      mlai.x{i} = mlai.x{i}(ii);
      mlai.y{i} = mlai.y{i}(jj);
      mlai.dat{i} = mlai.dat{i}(jj,ii);
    end
    dx = abs(median(diff(mlai.x{1})));
    dy = abs(median(diff(mlai.y{1})));
    mlai.als_lai = ones(size(mlai.dat{1}))*NaN;
    mlai.als_std = mlai.als_lai;
    for i = 1:length(mlai.x{1})
      for j = 1:length(mlai.y{1})
        ii = X >= mlai.x{1}(i)-dx/2 & X <= mlai.x{1}(i)+dx/2 & ...
             Y >= mlai.y{1}(j)-dy/2 & Y <= mlai.y{1}(j)+dy/2;
        tlai = LAI(ii);
        mlai.als_lai(j,i) =  nanmean(tlai(:));
        mlai.als_std(j,i) =  nanstd(tlai(:));
      end
    end
    eval(['save -v7.3 ~/Projects/Archived/BernJura/vegetation_res',num2str(res),'_cls.mat'])
end
eval(['load ~/Projects/Archived/BernJura/vegetation_res',num2str(res),'_cls.mat'])
clear lai slai hgt shgt
dlai = 0:0.25:7;
%MLAI = (MLAI+MLAI2)./2.*0.1;
for i = 1:length(dlai)-1;
  ii = MLAI >= dlai(i) & MLAI <= dlai(i+1) & VEH > 0;
  av_lai(i) = nanmean(LAI(ii));
  av_mlai(i) = nanmean(MLAI(ii));
  st_lai(i) = nanstd(LAI(ii))./sqrt(sum(ii(:)));
  st_mlai(i) = nanstd(MLAI(ii))./sqrt(sum(ii(:)));
end
errorbar(av_mlai,av_lai,st_lai);
xlabel('MODIS LAI');
ylabel('ALS LAI Proxy');
grid on
title('Comparison of LAI products : Berner Jura');
print -depsc2 /Users/morsdorf/Desktop/modis_lai_vs_mylai_both_onlyforest_onlyRTMAIN.eps
return

%jds = cal2jd(dat(:,3),dat(:,2),dat(:,1) + dat(:,4) / 24 + dat(:,5) / 60 / 24);
%jde = cal2jd(dat(:,3),dat(:,2),dat(:,1) + dat(:,6) / 24 + dat(:,7) / 60 / 24);
%[gpsweek,sows,rollover] = jd2gps(jds);
%[gpsweek,sowe,rollover] = jd2gps(jde);
%JDS = GPS*NaN;
dh = 100;
ht = [300:dh:1600];
%for i = 1:length(sows)
%    jj = find(GPS >= sows(i) & GPS <= sowe(i));
%    JDS(jj) = jds(i);
%    for j = 1:length(hth)-1
%        ii = find(GPS >= sows(i) & GPS <= sowe(i) & TER >= hth(j) & TER <= hth(j+1));
%        lai(j,i) = nanmean(LAI(ii));
%    end
%end
%figure;
%hlai1 = nanmean(lai(:,1:5)');
%hlai2 = nanmean(lai(:,6:8)');
%plot(ht,hlai1,'x-',ht,hlai2,'-*g')
%hold on
%plot(ht,lai(:,9),'o-r')

for j = 1:length(ht)
  ii = CLS == 1 & TER >= ht(j)-dh/2 & TER <= ht(j)+dh/2;
  slai(j) = nanstd(LAI(ii))./sqrt(sum(ii(:)));
  lai(j) = nanmean(LAI(ii));
  sfco(j) = nanstd(FCO(ii))./sqrt(sum(ii(:)));
  fco(j) = nanmean(FCO(ii));
  shgt(j) = nanstd(VEH(ii))./sqrt(sum(ii(:)));
  hgt(j) = nanmean(VEH(ii));
end
figure;
subplot(3,1,1)
errorbar(ht,lai,slai,'k');
grid on
xlabel('Altitude [m]');
ylabel('LAI Proxy');
title('Berner Jura - LAI proxy vs. Altitude, deciduous only'); 
%print -depsc2 ~/Projects/Archived/BernJura/laivsalti_decid.eps
subplot(3,1,2)
errorbar(ht,hgt,shgt,'k');
grid on
xlabel('Altitude [m]');
ylabel('Tree Height [m]');
title('Berner Jura - Tree height vs. Altitude, deciduous only'); 
subplot(3,1,3);
errorbar(ht,fco,sfco,'k');
grid on
xlabel('Altitude [m]');
ylabel('Fractional Cover [m]');
title('Berner Jura - Fractional Cover vs. Altitude, deciduous only'); 
print -depsc2 ~/Projects/Archived/BernJura/vegvsalti_decid.eps

for j = 1:length(ht)
  ii = CLS == 2 & TER >= ht(j)-dh/2 & TER <= ht(j)+dh/2;
  slai(j) = nanstd(LAI(ii))./sqrt(sum(ii(:)));
  lai(j) = nanmean(LAI(ii));
  sfco(j) = nanstd(FCO(ii))./sqrt(sum(ii(:)));
  fco(j) = nanmean(FCO(ii));
  shgt(j) = nanstd(VEH(ii))./sqrt(sum(ii(:)));
  hgt(j) = nanmean(VEH(ii));
end
figure;
subplot(3,1,1)
errorbar(ht,lai,slai,'k');
grid on
xlabel('Altitude [m]');
ylabel('LAI Proxy');
title('Berner Jura - LAI proxy vs. Altitude, conifers only'); 
%print -depsc2 ~/Projects/Archived/BernJura/laivsalti_decid.eps
subplot(3,1,2)
errorbar(ht,hgt,shgt,'k');
grid on
xlabel('Altitude [m]');
ylabel('Tree Height [m]');
title('Berner Jura - Tree height vs. Altitude, conifers only'); 
subplot(3,1,3);
errorbar(ht,fco,sfco,'k');
grid on
xlabel('Altitude [m]');
ylabel('Fractional Cover [m]');
title('Berner Jura - Fractional Cover vs. Altitude, conifers only'); 
print -depsc2 ~/Projects/Archived/BernJura/vegvsalti_conif.eps

return
load /Volumes/Data3/BernerJura/MyColorMap
%TER(isnan(TER)) = 9999;
%TER(TER==-9999) = NaN;
%TER = inpaint_nans(TER,4);
%TER(TER==9999) = NaN;
clear FCO VEH GPS
IMG1 = ones([size(LAI),3]);
IMG1(:,:,1) = JDS;
IMG1(:,:,2) = LAI;
IMG1(:,:,3) = TER;
for i = 1:3
    IMG1(:,:,i) = (IMG1(:,:,i)-min(min(IMG1(:,:,i))))/(max(max(IMG1(:,:,i))) - min(min(IMG1(:,:,i))));
end
return
strdat = {'FCT','FCO','LAI','VEH','TER','PDM','GPS'};
strt = {'Fractional Cover, all Angles.','Fractional Cover, small Angles.','Leaf area index','Average Canopy Height','Terrain','Point Density','GPSTime'};
shadmod(x(1:3:end),y(1:3:end),TER(1:3:end,1:3:end),LAI(1:1:end,1:1:end))
colormap(cmap)
swisstick
grid on
hcb = colorbar('South');
pos = get(hcb,'position');
set(hcb,'position',[pos(1)*3 pos(2)*0.95 pos(3)*0.6 pos(4)*0.7])
eval(['print -r900 -depsc2 ',strdat{i},'_res',num2str(res),'_BernerJura.eps'])
return
for i = 1:length(strdat)
    eval(['dat = ',strdat{i},';']);
    if i == 6
        dat = dat/9;
    end
    dat(TER < 400) = NaN;
    figure(1);clf
    imagesc(x,y,dat,'AlphaData',~isnan(dat));axis xy;axis equal;axis tight;
    orient landscape
    set(gca,'xtick',2560000:5000:2620000);
    swisstick;
    title(strt{i})
    if i == 1 | i == 2 | i == 3
        caxis([0 1])
        colormap(cmap)
    elseif i == 4
        colormap(jet)
        caxis([0 40])
    elseif i == 6
        colormap(jet)
        caxis([0 20]);
    else
        colormap(jet)
    end
    grid on
    hcb = colorbar('South');
    pos = get(hcb,'position');
    set(hcb,'position',[pos(1)*2.8 pos(2) pos(3)*0.6 pos(4)*0.7])
    set(gcf,'renderer','painters');
    eval(['print -r900 -depsc2 ',strdat{i},'_res',num2str(res),'_BernerJura.eps'])
end
