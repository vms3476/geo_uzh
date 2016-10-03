% script to compute biomass for swiss national park
load /Users/morsdorf/data/SNP/RASTER2010/dtmdsm2010.mat
chm = dsm;
chm.z = dsm.z - dtm.z;
%bio.dat = nlfilter(chm.z,[40 40],@biomass);
bio.dat = blockproc(chm.z,[40 40],@biomass);
bio.dat(bio.dat > 400) = NaN;
fun = @(block_struct) nanmean(block_struct.data(:));
bio.alt = blockproc(dtm.z,[40 40],fun);
res = 100;
alt = 1600:res:2500;
for i = 1:length(alt)
    ii = bio.alt >= alt(i)-res/2 & bio.alt <= alt(i)+res/2 & bio.dat > 10;
    mbio(i) = nanmean(bio.dat(ii));
    sbio(i) = nanstd(bio.dat(ii));
end
figure;
errorbar(alt,mbio,sbio);