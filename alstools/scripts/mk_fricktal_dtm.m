% script for placing all dtm's in a fricktal folder
cd /Users/morsdorf/data/KantonAargau/DTM-GeoTIFF
filez = dir('*.tif');

xi = [635,654];
yi = [248,270];

for i = 1:length(filez)
    nam = filez(i).name;
    xc = str2num(nam(5:7));
    yc = str2num(nam(12:14));
    if xc >= xi(1) & xc <= xi(2) & yc >= yi(1) & yc <= yi(2)
        disp(['Copying ',filez(i).name])
        unix(['cp ',filez(i).name,' Fricktal']);
    end
end
