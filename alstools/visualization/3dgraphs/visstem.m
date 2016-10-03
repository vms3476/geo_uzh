function [] = visstem(x,y,z,h,fname,col);
% function [] = visstem(x,y,z,h,fname,col);  
% function for stem visualization with osg
%   
% Felix Morsdorf, RSL Zurich, May 2005

  if nargin == 5
    col = [0.4 0.4 0];
  end
  ii = isnan(x);
  x(ii) = [];  y(ii) = [];  z(ii) = [];  h(ii) = [];
  
  figure
  plot3(x,y,z+h+2,'o','markersize',5,'color',col);
  hold on

  dx = 0.3;
  dy = 0.3;
  scl = 5;
  for i = 1:length(x)
    fill3([x(i)-dx,x(i)+dx,x(i)+(dx/scl),x(i)-(dx/scl),x(i)-dx], ...
          [y(i)-dy,y(i)+dy,y(i)+(dy/scl),y(i)-(dy/scl),y(i)-dy], ...
          [z(i),z(i),z(i)+h(i),z(i)+h(i),z(i)],col);
  end
  
  for i = 1:length(x)
    fill3([x(i)+dx,x(i)-dx,x(i)-(dx/scl),x(i)+(dx/scl),x(i)+dx], ...
          [y(i)-dy,y(i)+dy,y(i)+(dy/scl),y(i)-(dy/scl),y(i)-dy], ...
          [z(i),z(i),z(i)+h(i),z(i)+h(i),z(i)],col);
  end
  mat3d2osg(fname);
