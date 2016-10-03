function [tree] = readtlstree(fname)
% function [tree] = readtlstree(fname)
% loads tls trees from ASCII files, formats described in ISD of 3DVegLab project
% Input : filename
% Output : tree structure
  
% Felix Morsdorf, RSL 2012
  
fid = fopen(fname);
tree.base = str2num(fgetl(fid));
tree.bbox = str2num(fgetl(fid));
tree.rbas = str2num(fgetl(fid));
tree.numc = str2num(fgetl(fid));
for i = 1:tree.numc
  tree.cyl(i,:) = str2num(fgetl(fid));
end
fclose(fid);