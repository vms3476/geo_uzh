cc
if 0
    load dsm2010.mat;
    load dtm2010.mat;
    load snppoly;
    [X,Y] = meshgrid(dsm.x,dsm.y);
    [XI,YI] = meshgrid(1:length(dsm.x),1:length(dsm.y));
    ii = inpolygon(XI,YI,xi,yi);
    dsm.z(~ii) = 0;
    dsm.z = inpaint_nans(dsm.z,4);
    dsm.z(dsm.z==0)=NaN;
    dtm.z(~ii) = 0;
    dtm.z = inpaint_nans(dtm.z,4);
    dtm.z(dtm.z==0) = NaN;
    chm.z = dsm.z-dtm.z;
    dtm.x = dtm.x(1:end-1);
    dtm.z = dtm.z(:,1:end-1);
    dsm.y = dsm.y(2:end);
    dsm.z = dsm.z(2:end,:);
    save dtmdsm2010.mat
else
    load dtmdsm2010
    chm.z = dsm.z - dtm.z;
    clear dsm
    shadmod(dtm.x(1:1:end),dtm.y(1:1:end),dtm.z(1:1:end,1:1:end),chm.z(1:1:end,1:1:end));
    swisstick
    orient landscape
    wysiwyg
    colormap([[0.7,0.7,0.7];hsv(10)])
    caxis([0 30]);
    colorbar
    print -depsc2 -r700 /Users/morsdorf/Desktop/ForestGrowth/snp_chm2010.eps
end