function  [m,n] = dyad2(j, orientation)
% dyad2 -- Indices of entire j-th dyad, with spacified orientation
%                  of 2-d wavelet xform
%  Usage
%    [m,n] = dyad2(j, orientation);
%  Inputs
%    j     integer (between 0 and log2(N)-1 )
%    orientation  'h', 'v', 'd' standing for 
%                horizontal, verical and diagonal
%  Outputs
%    m,n    list of all 2-D indices of wavelet coeffs 
%
%  Needs functions dyad and dyadc
%
%  Example of use 
%
%  >> pict = MakeImage('StickFigure',128);
%  >> wf = MakeONFilter('Haar',1);
%  >> wpict = FWT2_PO(pict, 5, wf);
%  >> [diagx, diagy] = dyad2(6,'d');
%  >> diag_det = wpict(diagx, diagy);
%
   if strcmp( orientation, 'd')
        m = dyad(j);
        n = dyad(j);
   elseif strcmp(orientation, 'h')
        m = dyad(j);
        n = dyadc(j);
   elseif strcmp(orientation, 'v')
        m = dyadc(j);
        n = dyad(j);
   else
    	disp(sprintf('dyad2: I don''t recognize <<%s>>',orientation))
   		disp('Allowable Orientations are h, v, and d')
   end
%
% Copyright (c) 2001 Brani Vidakovic
%         
%   
% Addition to WaveLab Version 802
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    