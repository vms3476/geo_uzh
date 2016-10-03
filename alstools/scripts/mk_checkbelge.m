% script for testing coregistration of belgium ALS data with field data

load /Volumes/Data1/Belgium/kersselaerspleyn/kersselaerspleyn.mat
fds = shaperead(['/Volumes/Data1/Belgium/FieldData/' ...
                 'Hyperforest_KersselaerspleynKV+16Cirkels.shp']);
trees = importdata('/Volumes/Data1/Belgium/kersselaerspleyn/Products/trees.csv',',',1); 
myimage(chm);
hold on
ii = [fds(:).DBH_mm] > 70;
x = [fds(ii).X];y = [fds(ii).Y];
hp(1) = plot(x,y,'.k')
hold on
hp(2) = plot(trees.data(:,1),trees.data(:,2),'.','color',[0.7 0.7 0.7]);
axis([min(x) max(x) min(y) max(y)])
mytitle(['Tree locations from ALS (gray) and field measurements (black), DBH > 70 ' ...
         'mm']);
set(hp, 'markersize', 3);
print -depsc2 /Users/morsdorf/Desktop/treeloc_kers.eps