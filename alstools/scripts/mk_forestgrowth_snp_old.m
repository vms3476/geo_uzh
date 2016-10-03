% script for producing a map of the forest growth in switzerland
%load old modelss
cc
if 0
cd /Users/morsdorf/Desktop/ForestGrowth
load dsm;dsm02.z = dsm02.z';
load dtm;dtm02.z = dtm02.z';

load dsm2010
load dtmf2010

disp('Load completed...')
ras = dtm02;
[ras.X,ras.Y] = meshgrid(ras.x,ras.y);
ras.dtm02 = ras.z;
[dsm02.X,dsm02.Y] = meshgrid(dsm02.x,dsm02.y);
ras.dsm02 = interp2(dsm02.x,dsm02.y,dsm02.z,ras.X,ras.Y,'cubic');
disp('1st interp completed...')
clear dsm02 dtm02
ii = isnan(ras.dsm02);
ras.dtm02(ii) = -9999;
ras.dtm02 = inpaint_nans(ras.dtm02,4);
ras.dtm02(ii) = NaN;
disp('Inpaint completed...')
ii = dtm.x >= min(ras.x) & dtm.x <= max(ras.x);
jj = dtm.y >= min(ras.y) & dtm.y <= max(ras.y);
dtm.x = double(dtm.x(ii));dtm.y = double(dtm.y(jj));dtm.z = double(dtm.z(jj,ii));
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
dtm.z = inpaint_nans(dtm.z,4);
ras.dtm10 = interp2(dtm.X,dtm.Y+1,double(dtm.z),ras.X,ras.Y,'cubic');
disp('3rd interp completed...')


ii = dsm.x >= min(ras.x) & dsm.x <= max(ras.x);
jj = dsm.y >= min(ras.y) & dsm.y <= max(ras.y);
dsm.x = double(dsm.x(ii));dsm.y = double(dsm.y(jj));dsm.z = double(dsm.z(jj,ii));
dsm.z = inpaint_nans(dsm.z,4);
disp('Inpaint completed...')

[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
ras.dsm10 = interp2(dsm.X,dsm.Y+1,dsm.z,ras.X,ras.Y,'cubic');
disp('2nd interp completed...')
s = sysshift(ras.dtm02(750:1250,2000:3000),ras.dtm10(750:1250,2000:3000));
save ras.mat 
end
load ras.mat
ras.dtm10 = interp2(dtm.X,dtm.Y-0.5,dtm.z,ras.X,ras.Y,'cubic');
ras.dsm10 = interp2(dsm.X,dsm.Y-0.5,dsm.z,ras.X,ras.Y,'cubic');
s = sysshift(ras.dtm02(750:1250,2000:3000),ras.dtm10(750:1250,2000:3000));
ras.dsm10(isnan(ras.dsm02)) = NaN;
ras.dtm10(isnan(ras.dsm02)) = NaN;
ras.chm02 = ras.dsm02-ras.dtm02;
ras.chm10 = ras.dsm10-ras.dtm10;
for i = 1
    clf
    if i==1
        pl_dtm = 1;
    else
        pl_dtm = 0;
    end
    if pl_dtm
        ras.diff = ras.dtm10-ras.dtm02;
    else
        ras.diff = ras.chm10-ras.chm02;
        h = fspecial('average',5);
        ras.diff = imfilter(ras.diff,h);
    end
    load polygon2
    ypol1(1) = ypol1(1)-20;
    ypol2(2) = ypol2(2)-20;
    ii = inpolygon(ras.X(:),ras.Y(:),xpol1([1:3,1]),ypol1([1:3,1]));
    ras.diff(ii) = NaN;
    ii = inpolygon(ras.X(:),ras.Y(:),xpol2([1:3,1]),ypol2([1:3,1]));
    ras.diff(ii) = NaN;
    ras.diff = interp2(ras.X,ras.Y,ras.diff,dtm.X,dtm.Y);
    ii = isnan(ras.diff);
    dtm.z(ii) = NaN;clear dsm;
    %NDTM.x = linspace(dtm.x(1),dtm.x(end),length(dtm.x)*1.5);
    %NDTM.y = linspace(dtm.y(1),dtm.y(end),length(dtm.y)*1.5);
    %keyboard
    %[NDTM.X,NDTM.Y] = meshgrid(NDTM.x,NDTM.y);
    %NDTM.Z = interp2(dtm.X,dtm.Y,dtm.z,NDTM.X,NDTM.Y,'cubic');
    %NDTM.diff = interp2(dtm.X,dtm.Y,ras.diff,NDTM.X,NDTM.Y,'cubic');
    %clear dtm ras
    %shadmod(NDTM.x,NDTM.y,NDTM.Z,NDTM.diff);
    shadmod(dtm.x,dtm.y,dtm.z,ras.diff);
    swisstick
    colormap(flipud(bipolar(256,0.7)))
    if pl_dtm
        caxis([-2 2])
    else
        caxis([-3 3])
    end
    cax = caxis;
    %brighten(0.5)
    view(-22.62,90)
    hl = findobj('type','light');
    [az,el] = lightangle(hl);
    lightangle(hl,az+22.62,el);
    %set(gca,'visible','off');
    %set(hcb,'rotation',-22.62);
    orient landscape
    set(gcf,'PaperType','A0');
    set(gca,'position',[-0.06 -0.11 0.77*1.5 0.77*1.5])
    title('Ofenpass Vegetation Change 2002 - 2010')
    hcb = colorbar('NorthOutside');
    set(hcb,'position',[0.075 0.69 0.17 0.015],'xaxislocation','bottom','xtick',[min(cax):0.5:max(cax)],'fontsize',6,'box','on')
    set(gcf,'renderer','zbuffer')
    if pl_dtm
        print -dtiff -r1200 /Users/morsdorf/Desktop/ForestGrowth/diff_dtm.tif
    else
        print -dtiff -r1200 /Users/morsdorf/Desktop/ForestGrowth/diff_chm.tif
    end
end