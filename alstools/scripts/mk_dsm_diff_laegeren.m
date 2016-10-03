% script for producing dsm differences for laegeren
cc
%cd /Volumes/lidarlab/data/Laegeren/Spring_2010/DSM_First_Echo
%filez = dir('*.coo');
%dsm10 = xyz2ras(filez);
%save DSM_2010_Leafoff dsm10

%load /Volumes/lidarlab/data/Laegeren/Summer_2010/dsm_fe/dsm.mat

load /Volumes/lidarlab/data/Laegeren/Spring_2010/DSM_First_Echo/DSM_2010_Leafoff.mat
dsm = dsm10;
if 1
    are = [min(dsm.x) max(dsm.x) min(dsm.y) max(dsm.y)];
    
    raw = getrawlas('/Volumes/lidarlab/data/Laegeren/2014_ZH/03_LAS_LV03_LN02/',are);
    ii = raw.rnnr == 11 | raw.rnnr == 21 | raw.rnnr == 31 | raw.rnnr == 41 | raw.rnnr == 51 | raw.rnnr == 61; 
    
    rawf = subsetraw(raw,ii);
    save /Users/morsdorf/Desktop/raw_laegeren_leafoff2014 rawf
    
    dsm14 = raw2ras([rawf.x,rawf.y,rawf.z],0.5,0.5,'dsm');
    save /Users/morsdorf/Desktop/dsm_laegeren_leafoff2014 dsm14
else
    load /Users/morsdorf/Desktop/dsm_laegeren_leafoff2014
end
 ii = dsm.x >= min(dsm14.x) & dsm.x <= max(dsm14.x);
 jj = dsm.y >= min(dsm14.y) & dsm.y <= max(dsm14.y);
 dsm.x = dsm.x(ii);
 dsm.y = dsm.y(jj);
 dsm.z = dsm.z(jj,ii);

fun = @(x) nanmean(x(:));
dsm.z = nlfilter(dsm.z,[6 6],fun);
dsm14.z = nlfilter(dsm14.z,[6 6],fun);
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
[dsm14.X,dsm14.Y] = meshgrid(dsm14.x,dsm14.y);
dz = dsm14.z - interp2(dsm.X,dsm.Y,dsm.z,dsm14.X,dsm14.Y);
ddsm = dsm14;
ddsm.z = dz;
clf
myimage(ddsm)
caxis([-5 5]);
colorbar
title('DSM Height Change 2010 - 2014, leaf-on');
