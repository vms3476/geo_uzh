% convert SHP to mat
cd /Users/morsdorf/Desktop/Shapefiles/
filez = dir('*.shp');
for i = 1:length(filez)
    disp([filez(i).name(1:end-4), ' = shaperead(filez(i).name);']);
    eval([filez(i).name(1:end-4), ' = shaperead(filez(i).name);']);
end
