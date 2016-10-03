% script for making 3d vis of lidar data from Kt. Aargau
if 0
cd /Users/morsdorf/Desktop/Thierstein3d

raw_on = readlas('640000_260000_on.las');

raw_off = readlas('640000_260000_off.las');

raw_on.x = raw_on.x - 640000;
raw_on.y = raw_on.y - 260000;

raw_off.x = raw_off.x - 640000;
raw_off.y = raw_off.y - 260000;


hl = pldiscret(raw_on,raw_off,'none');

mat3d2osg('Thierstein_onoff')
close;
figure;
end
raw_on.int(raw_on.int > 100) = 100;
raw_off.int(raw_off.int > 100) = 100;
myscatter3(raw_on.x,raw_on.y,raw_on.z,raw_on.int,gray(256),3);

myscatter3(raw_off.x,raw_off.y,raw_off.z,raw_off.int,gray(256),3);