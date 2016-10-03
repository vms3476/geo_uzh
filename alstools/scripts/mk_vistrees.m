% script for visualizing the trees from Laegeren
if 1
  cc
  cd /Users/morsdorf/Projects/Active/3DVegLab/ISD/Laegeren
  load tree_matched.mat
  cd trees_TLS
  filez = dir('*.asc');
  for i = 1:length(filez);
    tree{i} = readtlstree(filez(i).name);
    tree{i}.base(1) = tree{i}.base(1) + 600000;
    tree{i}.base(2) = tree{i}.base(2) + 200000;
    BASE(i,:) = tree{i}.base;
    if tree{i}.base(1) < 669770
      figure(1);
      %hp1(i) = plottlstrees(tree{i});
      axis equal
      view(3)
      swisstick
    else
      figure(2);
      %hp2(i) = plottlstrees(tree{i});
      axis equal
      view(3)
      swisstick
    end
  end
end

figure(1);
%ax = axis;
ax = [669654 669709 259045 259092];
ii = tree_final.x >= ax(1) & tree_final.x <= ax(2) &  ...
     tree_final.y >= ax(3) & tree_final.y <= ax(4) ;
trees1 = subsetraw(tree_final,ii);
jj = [];
for i = 1:length(trees1.x)
  ii = BASE(:,1) >= trees1.x(i) - trees1.dia_WE(i)/1.8 & ...
       BASE(:,1) <= trees1.x(i) + trees1.dia_WE(i)/1.8 & ...
       BASE(:,2) >= trees1.y(i) - trees1.dia_NS(i)/1.8 & ...
       BASE(:,2) <= trees1.y(i) + trees1.dia_NS(i)/1.8;
  if sum(ii) ~= 0
    jj = [jj,i];
  end
end
trees1 = subsetraw(trees1,jj);     
ht1 = treemodel3(trees1.x,trees1.y,trees1.z,trees1.treeheight,trees1.crown_length, ...
                 [trees1.dia_WE,trees1.dia_NS],ones(size(trees1.x))*2);
camlight
set(gca,'zlim',[min(trees1.z) max(trees1.z)+max(trees1.treeheight)]);
% clean up
print -depsc2 -r600 plot1_als3d.eps
print -dpng -r600 plot1_als3d.png
view(2)
print -dpng -r600 plot1_als2d.png
print -depsc2 -r600 plot1_als2d.eps
figure(2)
%ax = axis;
ax = [669818 669867 259010 259067];
ii = tree_final.x >= ax(1) & tree_final.x <= ax(2) &  ...
     tree_final.y >= ax(3) & tree_final.y <= ax(4) ;
trees2 = subsetraw(tree_final,ii); 
jj = [];
for i = 1:length(trees2.x)
  ii = BASE(:,1) >= trees2.x(i) - trees2.dia_WE(i)/1.8 & ...
       BASE(:,1) <= trees2.x(i) + trees2.dia_WE(i)/1.8 & ...
       BASE(:,2) >= trees2.y(i) - trees2.dia_NS(i)/1.8 & ...
       BASE(:,2) <= trees2.y(i) + trees2.dia_NS(i)/1.8;
  if sum(ii) ~= 0
    jj = [jj,i];
  end
end
trees2 = subsetraw(trees2,jj);     
ht2 = treemodel3(trees2.x,trees2.y,trees2.z,trees2.treeheight,trees2.crown_length, ...
                 [trees2.dia_WE,trees2.dia_NS],ones(size(trees2.x))*2);
camlight
set(gca,'zlim',[min(trees2.z) max(trees2.z)+max(trees2.treeheight)]);
print -depsc2 -r600 plot2_als3d.eps
print -dpng -r600 plot2_als3d.png
view(2)
print -dpng -r600 plot2_als2d.png
print -depsc2 -r600 plot2_als2d.eps