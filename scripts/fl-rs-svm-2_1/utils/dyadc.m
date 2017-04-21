function i = dyadc(j)
% dyadc -- Indices preceeding the entire j-th dyad of 1-d wavelet xform
%  Usage
%    i = dyadc(j);
%  Inputs
%    j     integer
%  Outputs
%    i    list of all indices of wavelet coeffts starting at 1
%         and complemenying j-th level
%
    i = 1:2^(j);
%
% Copyright (c) 2001 Brani Vidakovic
%     
    
    
%   
% Addition to WaveLab Version 802
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
