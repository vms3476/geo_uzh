% script to compute distances between echo types for Laegeren
cd /Users/morsdorf/Desktop
raw = readlas('/Users/morsdorf/data/Laegeren/2010/RAW_SUMMER/data.las');

rnnr = double(raw.rnnr);

dz = diff(raw.z);
drnnr = diff(rnnr);
unique(drnnr);
dz(dz < -10) = -10;
ii = drnnr == 1 & dz < 0;
clf;myscatter3(raw.x(ii),raw.y(ii),raw.z(ii),dz(ii),hsv(256),2);
mat3d2osg('flx_twr_2010_summer_echodist.osg');



raw = readlas('/Users/morsdorf/data/Laegeren/2010/RAW_SPRING/data.las');

rnnr = double(raw.rnnr);

dz = diff(raw.z);
drnnr = diff(rnnr);
unique(drnnr);

dz(dz < -10) = -10;
ii = drnnr == 1 & dz < 0;
clf;myscatter3(raw.x(ii),raw.y(ii),raw.z(ii),dz(ii),hsv(256),2);
mat3d2osg('flx_twr_2010_spring_echodist.osg');

