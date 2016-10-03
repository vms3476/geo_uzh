% script to visualize/convert the stem models of Tharandt
cd /Users/morsdorf/Projects/Active/3DVegLab/3D_RekTLS/Tree_Models_Wavefront_objects
filez = dir('*.obj');

%for i = 1:length(filez)
%    tree{i} = read_wobj(filez(i).name);
%end
%save trees.mat tree
load trees