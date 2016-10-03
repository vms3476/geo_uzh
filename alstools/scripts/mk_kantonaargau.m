cc
cd /Users/morsdorf/Projects/Active/KantonAargau/Daten/Einzelbaumerkennung/Polygons
filez = dir('6*.shp');
cmap = hsv(64);
for i = 1:length(filez)
    if i ~= 4
        S = shaperead(filez(i).name);
    else
        load 669000_259000_poly.mat
    end
        colind = fix(rand(length(S))*64)+1;
    for j = 1:length(S)
        plot(S(j).X,S(j).Y,'-','color',cmap(colind(j),:));
        hold on
    end
end
cd /Users/morsdorf/Projects/Active/KantonAargau/Daten/Einzelbaumerkennung/Points
filez = dir('6*.shp');
for i = 1:length(filez)
    S = shaperead(filez(i).name);
    for j = 1:length(S)
        plot(S(j).X,S(j).Y,'.k');
        hold on
    end
end