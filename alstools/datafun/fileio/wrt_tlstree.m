function [] = wrt_tlstree(tree,fname);
% function [] = wrt_tlstree(tree,fname);
% writes TLS based tree reconstruction (stem and branch cylinders) to ASCII file
% format is described in ISD document of the 3DVegLab project
% tree is a structure derived from reading .obj file with read_wobj (from filecentral@mathworks).

% Felix Morsdorf, March 2013

% open ASCII file for writing
fid = fopen(fname,'w');

% determine base point

