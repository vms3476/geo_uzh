% script for loading and displaying the data from neuchatel
close all
cd /Users/morsdorf/data/1164313
if 0
    dtm = asciigrid2mat('1164313.asc');
    raw = readlas('1164313.las');
    [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
    raw.tz = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
    dtm.oz = dtm.z;
    [m,n] = size(dtm.z);
    dtm.z = imresize(dtm.z,[m/2,n/2]);
    dtm.ox = dtm.x;
    dtm.oy = dtm.y;
    dtm.x = dtm.x(1:2:end);
    dtm.y = dtm.y(1:2:end);
    dsm = raw2ras([raw.x,raw.y,raw.z],dtm,1,'dsm');
    cls = raw2cls([raw.x,raw.y,double(raw.Classification)],dtm,1);
    chm = dtm;
    chm.z = dsm.z-dtm.z;
    ii = chm.x >= min(raw.x) & chm.x <= max(raw.x);
    jj = chm.y >= min(raw.y) & chm.y <= max(raw.y);
    chm.z = chm.z(jj,ii);
    chm.x = chm.x(ii);
    chm.y = chm.y(jj);
    ii = cls.x >= min(raw.x) & cls.x <= max(raw.x);
    jj = cls.y >= min(raw.y) & cls.y <= max(raw.y);
    cls.z = cls.z(jj,ii);
    cls.x = cls.x(ii);
    cls.y = cls.y(jj);
    chm.z = inpaint_nans(chm.z,4);
    save -v7.3 data
else
    %load data
end
figure(1);clf
myimage(chm)
ii = locmax(chm.z,9,0.7,3);
hold on
jj = cls.z == 6 | cls.z == 2 | cls.z > 6;
ii(jj(:)) = 0;
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);
plot(chm.X(ii),chm.Y(ii),'.k','markersize',3)
loc = [chm.X(ii),chm.Y(ii),chm.z(ii)];
raw.oz = raw.z;
raw.z = raw.tz;
trees = segtree(raw,loc,3,4);