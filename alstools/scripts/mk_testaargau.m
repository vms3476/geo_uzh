% script for testing data generated with different thresholds
if 0
dtm3 = xyz2ras('/Users/morsdorf/data/TEST_KTAARGAU/DTM/SW3/32413000_5232000.asc');

dtm5 = xyz2ras('/Users/morsdorf/data/TEST_KTAARGAU/DTM/SW5/32413000_5232000.asc');


dsm3 = xyz2ras('/Users/morsdorf/data/TEST_KTAARGAU/DSM/SW3/32413000_5232000.asc');

dsm5 = xyz2ras('/Users/morsdorf/data/TEST_KTAARGAU/DSM/SW5/32413000_5232000.asc');
end

str = {'DTM','DSM'};
for i = 1:2
    figure(i);clf
    if i == 1
        diff = dtm3.z-dtm5.z;
    else
        diff = dsm3.z-dsm5.z;
    end
    myimage(dtm3.x,dtm3.y,diff);
    caxis([-0.5 0.5])
    colormap(ocean2);
    colorbar
    swisstick;
    title([str{i},' Difference 3 - 5']);
    print('-depsc2',[str{i},'_diff.eps']);
end
figure(3);clf
subplot(1,2,1)
shadmod(dtm3)
title('DTM 3');
subplot(1,2,2);
shadmod(dtm5)
title('DTM 5');
colormap(gray);
%raw5 = readlas('/Users/morsdorf/data/TEST_KTAARGAU/Las/SW5/32413000_5232000.las');

%raw3 = readlas('/Users/morsdorf/data/TEST_KTAARGAU/Las/SW3/32413000_5232000.las');
