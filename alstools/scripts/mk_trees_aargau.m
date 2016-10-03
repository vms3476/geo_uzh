cc;
[hdr,data] = hdrload('/Users/morsdorf/Desktop/KantonAargauSingleTrees/trees.csv');
[hdr2,data2] = hdrload('/Users/morsdorf/Desktop/KantonAargauSingleTrees/trees_ext.csv');
for i = 1:length(data(:,1))
    ii = find(data(i,1) == data2(:,1));
    data(i,4:6) = data2(ii,[3 8 12]);
end
perc = 10;
ii = data(:,4) > perc;
hp(1) = scatter(data(ii,2),data(ii,3),round(data(ii,6)*3),'r');
hold on
hp(2) = scatter(data(~ii,2),data(~ii,3),round(data(~ii,6)*3),'g');
hold on
legend(hp,{['> ',num2str(perc),'% conifers'],['< ',num2str(perc),'% conifers']});
%ht = text(data(:,2),data(:,3),num2str([data(:,1)]));
axis([0 100 0 100])
plot([0,100],[0,100],'-k','linewidth',2);
xlabel('Field measured')
ylabel('LiDAR measured')
title('Number of trees per stand unit')
ii = abs(data(:,2)-data(:,3)) < 40;
print -depsc2 tree_kantonaargau.eps
%add_regress(data(ii,2),data(ii,3));
%print -depsc2 tree_kantonaargau_regress.eps