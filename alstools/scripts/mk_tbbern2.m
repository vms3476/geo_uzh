% script for producing graphics for the technical report in Bern
cc
cd /Volumes/Data3/BernerJura/Punktwolke_klassiert/
% make global histograms of shots and echos
load /Volumes/Data3/BernerJura/stats.mat


% make histograms of point density
figure(1)
clf
histshots = ones(length(stats),30);
histechos = ones(length(stats),30);
histgnd = ones(length(stats),30);
histncl = ones(length(stats),30);
histovr = ones(length(stats),30);
for i = 1:length(stats)
    histechos(i,:) = stats(i).histechos;
    histshots(i,:) = stats(i).histshots;  
    histgnd(i,:) =   stats(i).histclass(2,:);  
    histncl(i,:) =   stats(i).histclass(1,:);  
    histovr(i,:) =   stats(i).histclass(12,:);  
end

subplot(2,1,1)
N = sum(histshots)./sum(sum(histshots))*100;
X = 1:30;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 20]);
xlabel('Schuesse pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = sum(histshots);
title(['Histogramm Schussdichte alle Kacheln (n = ',num2str(length(stats)), ...
       ') , Mittel: ',num2str(mean([stats.meanshots]))]);

subplot(2,1,2)
N = sum(histechos)./sum(sum(histechos))*100;
X = 1:30;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 20]);
xlabel('Echos pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = sum(histechos);
title(['Histogramm Echodichte alle Kacheln (n = ',num2str(length(stats)), ...
       ') , Mittel: ',num2str(mean([stats.meanechos]))]);
print -r600 -depsc2 PunktdichteHistogramAlleKacheln.eps


figure(2)
N = sum(histgnd)./sum(sum(histgnd))*100;
X = 1:30;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 55]);
xlabel('Bodenpunkte pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
title(['Histogramm Bodenpunkte alle Kacheln, n = ',num2str(length(stats))]);
print -r600 -depsc2 BodenpunkteHistogramAlleKacheln.eps

figure(3)
N = sum(histovr)./sum(sum(histovr))*100;
X = 1:30;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 75]);
xlabel('Ueberlappungspunkte pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
title(['Histogramm Ueberlappungspunkte alle Kacheln, n = ',num2str(length(stats))]);
print -r600 -depsc2 OverlapHistogramAlleKacheln.eps


figure(4)
N = sum(histncl)./sum(sum(histncl))*100;
X = 1:30;
bar(X,N);
set(gca,'xtick',[0:2:30],'ylim',[0 75]);
xlabel('Ueberlappungspunkte pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
title(['Histogramm unklassierte Punkte alle Kacheln, n = ',num2str(length(stats))]);
print -r600 -depsc2 NotClassifiedHistogramAlleKacheln.eps

% make map of point density
if 0
    cd /Volumes/Data3/BernerJura/Punktwolke_klassiert/
    filez = dir('/Volumes/Data3/BernerJura/Punktwolke_klassiert/*_*.mat');
    res = 3;
    x = 2555000:res:2621000;
    y = 1211000:res:1245000;
    [X,Y] = meshgrid(x,y);
    NUMSHOTS = ones(size(X))*NaN;
    NUMECHOS = NUMSHOTS;
    %NUMECHOV = NUMSHOTS;
    NUMGND = NUMSHOTS;
    NUMVEG = NUMSHOTS;
    NUMBLD = NUMSHOTS;
    NUMOVR = NUMSHOTS;
    %NUMUCL = NUMSHOTS;
    
    for i = 1:length(filez)
        disp(filez(i).name);
        load(filez(i).name);
        %numshots(isnan(numshots)) = 0;
        %numechos(isnan(numechos)) = 0;
        %numechov(isnan(numechov)) = 0;
        %numclass(isnan(numclass)) = 0;        
        %xg = min(xg)+(res/2):res:max(xg)-(res/2);
        %yg = min(yg)+(res/2):res:max(yg)-(res/2);
        [XG,YG] = meshgrid(xg,yg);

        if min(size(XG)) > 2
            II = find(x >= min(xg)  &  x <= max(xg));
            JJ = find(y >= min(yg)  &  y <= max(yg));
            fun = @(block_struct) nanmean(block_struct.data);
            numshots = colfilt(numshots,[3 3],'sliding',@nanmean);
            numechos = colfilt(numechos,[3 3],'sliding',@nanmean);
            %numechov = colfilt(numechov,[3 3],'sliding',@nanmean);
            numclass2 = numclass;
            for o = [2,3,5,6,12]
                numclass2(:,:,o) = colfilt(numclass(:,:,o),[3 3],'sliding',@nanmean);
            end
            numclass = numclass2;clear numclass2;
            if res > 5
                for k = 1:length(II)
                    for l = 1:length(JJ)
                        ii = find(xg >= x(II(k))-(res/2)  &  xg <= x(II(k))+(res/2));
                        jj = find(yg >= y(JJ(l))-(res/2)  &  yg <= y(JJ(l))+(res/2));
                        NUMSHOTS(JJ(l),II(k)) = nanmean(nanmean(numshots(jj,ii)));
                        NUMECHOS(JJ(l),II(k)) = nanmean(nanmean(numechos(jj,ii)));
                        NUMGND(JJ(l),II(k)) = nanmean(nanmean(nanmean(numclass(jj,ii,2))));
                        NUMVEG(JJ(l),II(k)) = nanmean(nanmean(nanmean(numclass(jj,ii,[3,5]))));
                        NUMBLD(JJ(l),II(k)) = nanmean(nanmean(nanmean(numclass(jj,ii,6))));
                        %NUMUCL(JJ(l),II(k)) = nanmean(nanmean(nanmean(numclass(jj,ii,1))));
                        NUMOVR(JJ(l),II(k)) = nanmean(nanmean(nanmean(numclass(jj,ii,12))));
                    end
                end
            else
                if min(size(XG)) > 2;
                    numveg = numclass(:,:,3)+numclass(:,:,5);
                    NUMSHOTS(JJ,II) = interp2(XG,YG,numshots,X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    NUMECHOS(JJ,II) = interp2(XG,YG,numechos,X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    %NUMECHOV(JJ,II) = interp2(XG,YG,numechov,X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    NUMGND(JJ,II) = interp2(XG,YG,numclass(:,:,2),X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    NUMVEG(JJ,II) = interp2(XG,YG,numveg,X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    NUMBLD(JJ,II) = interp2(XG,YG,numclass(:,:,6),X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    %NUMUCL(JJ,II) = interp2(XG,YG,numclass(:,:,1),X(JJ,II),Y(JJ,II),'*nearest',-1000);
                    NUMOVR(JJ,II) = interp2(XG,YG,numclass(:,:,12),X(JJ,II),Y(JJ,II),'*nearest',-1000);
                end
            end
        end
    end
    save maps.mat NUMOVR NUMSHOTS NUMECHOS NUMGND NUMBLD NUMVEG x y
end
load maps 
mask = NUMGND > 0 | NUMBLD > 0 | NUMVEG > 0 | NUMOVR > 0;
load ../cmap
figure(1);clf
subplot(2,1,1)
X = 0:30;
nshots = NUMSHOTS(mask);
[N,BIN] = histc(nshots(:),X);
N = N./sum(N);
bar(X,N*100,'histc');
set(gca,'xtick',[0:1:30],'ylim',[0 20]);
xlabel('Schuesse pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = sum(NUMSHOTS(:));
title(['Histogramm Schussdichte alle Kacheln (n = ',num2str(length(stats)), ...
       ') , Mittel: ',num2str(mean([nshots(:)]))]);

subplot(2,1,2)

X = 0:30;
nechos = NUMECHOS(mask);
[N,BIN] = histc(nechos(:),X);
N = N./sum(N);
bar(X,N*100,'histc');
set(gca,'xtick',[0:1:30],'ylim',[0 20]);
xlabel('Echos pro Quadratmeter');
ylabel('Anzahl Pixel [%]');
NH = sum(NUMECHOS(:));
title(['Histogramm Echodichte alle Kacheln (n = ',num2str(length(stats)), ...
       ') , Mittel: ',num2str(mean([nechos(:)]))]);
%print -r600 -depsc2 PunktdichteHistogramAlleKacheln3m.eps

figure(2);clf
imagesc(x,y,NUMSHOTS,'AlphaData',[~isnan(NUMSHOTS) & NUMECHOS > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000);
swisstick;
title('Schussdichte Berner Jura')
caxis([0 20])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 18],'xtick',[0:2:18],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Schussdichte [m^{-2}]');

%print -r600 -depsc2 SchussdichteBernerJura.eps
%wrt_geotiff_CH('Schussdichte',x,y,NUMSHOTS);
figure(3);clf
imagesc(x,y,NUMECHOS,'AlphaData',[~isnan(NUMECHOS) & NUMECHOS > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Echodichte Berner Jura')
caxis([0 20])
hax = gca;
grid on
colormap([cmap])
hcb(2) = colorbar('South');
set(hcb,'xlim',[0 18],'xtick',[0:2:18],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb(2));xlabel('Echodichte [m^{-2}]');
print -r600 -depsc2 EchodichteBernerJura.eps
axes(hax);hold on
load /Volumes/Data3/BernerJura/buildings.mat
plot(MEX,MEY,'ok','markersize',10);
plot(MEX,MEY,'xw','markersize',10);
str = {'Bevilard','Cortebert','Moutier','Orvin','Vauffelin'};
text(MEX+500,MEY,str,'fontsize',12,'fontweight','bold');
print -r600 -depsc2 EchodichteBernerJura_Villages.eps
%wrt_geotiff_CH('Echodichte',x,y,NUMECHOS);

figure(4);clf
imagesc(x,y,NUMECHOS,'AlphaData',[~isnan(NUMECHOS) & NUMECHOS < 4 & NUMECHOS > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Gebiete mit weniger als 4 Echos/m^2')
caxis([0 4])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 4],'xtick',[0:0.5:4],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Echodichte [m^{-2}]');
%print -r600 -depsc2 Weniger4EchosBernerJura.eps





figure(5);clf
imagesc(x,y,NUMSHOTS,'AlphaData',[~isnan(NUMSHOTS) & NUMSHOTS < 4 & NUMSHOTS > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Gebiete mit weniger als 4 Schuesse/m^2')
caxis([0 4])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 4],'xtick',[0:0.5:4],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Schussdichte [m^{-2}]');
%print -r600 -depsc2 Weniger4SchussBernerJura.eps


figure(6);clf
imagesc(x,y,NUMGND,'AlphaData',[~isnan(NUMGND) & NUMGND > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Bodenpunkte pro Quadratmeter')
caxis([0 10])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 10],'xtick',[0:1:10],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Bodenpunkte [m^{-2}]');
%print -r600 -depsc2 BodenpunkteBernerJura.eps
%wrt_geotiff_CH('Bodenpunkte',x,y,NUMGND);
figure(7);clf
imagesc(x,y,NUMVEG,'AlphaData',[~isnan(NUMVEG) & NUMVEG > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Vegetationspunkte pro Quadratmeter')
caxis([0 10])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 10],'xtick',[0:1:10],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Vegetationspunkte [m^{-2}]');
%print -r600 -depsc2 VegetationspunkteBernerJura.eps
%wrt_geotiff_CH('Vegetationspunkte',x,y,NUMVEG);

figure(8);clf
imagesc(x,y,NUMBLD,'AlphaData',[~isnan(NUMBLD) & NUMBLD > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Gebaeudepunkte pro Quadratmeter')
caxis([0 10])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 10],'xtick',[0:1:10],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Gebaeudepunkte [m^{-2}]');
print -r600 -depsc2 GebaeudepunkteBernerJura.eps
wrt_geotiff_CH('Gebaeudepunkte',x,y,NUMBLD);
return
figure(9);clf
imagesc(x,y,NUMOVR,'AlphaData',[~isnan(NUMOVR) & NUMOVR > 0]);axis xy;axis equal;axis tight;
orient landscape
set(gca,'xtick',2560000:5000:2620000)
swisstick;
title('Ueberlappungspunkte pro Quadratmeter')
caxis([0 10])
grid on
hcb = colorbar('South');
colormap([cmap])
pos = get(hcb,'position');
set(hcb,'xlim',[0 10],'xtick',[0:1:10],'position',[pos(1)+0.4 pos(2) pos(3)-0.4 pos(4)]);
axes(hcb);xlabel('Ueberlappungspunkte [m^{-2}]');
print -r600 -depsc2 UeberlappungspunkteBernerJura.eps
%wrt_geotiff_CH('Ueberlappungspunkte',x,y,NUMOVR);
