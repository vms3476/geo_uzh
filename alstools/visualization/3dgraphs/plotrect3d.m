function [h,pos] = plotrect3d(are,dtm);
% function [h,pos] = plotrect3d(are,dtm);
% plots a rectangle defined by four coordinates as interpolated line on DTM or DSM model
% or delivers coordinates of line segements, if used with two output arguments
% Inputs :
%          are = [xmin xmax ymin ymax];
%          dtm = structure containing dtm.x,dtm.y,dtm.z (optional as well dtm.X and dtm.y as 
%                obtained by meshgrid
% Outputs:
%          h   = handle to 3d lineobject
%          pos = structure containing line coordinates (pos.x,pos.y,pos.z)

% Felix Morsdorf, RSL Zuerich, May 2013

  % compute spacing
  dx = median(diff(dtm.x));
  dy = median(diff(dtm.y));
  
  % get length of segments
  xd = fix(are(2)-are(1)/dx);
  yd = fix(are(4)-are(3)/dy);
  
  % define 2d polygon
  if min(size(are)) == 1 % bbox
      x = [ones(1,yd)*are(1) linspace(are(1),are(2),xd) ones(1,yd)*are(2) linspace(are(2),are(1),xd)];
      y = [linspace(are(3),are(4),yd) ones(1,xd)*are(4) linspace(are(4),are(3),yd)  ones(1,xd)*are(3)];
  end
  
  % interpolate z coordinate
  if ~isfield(dtm,'X')
      [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
  end
  
  z = interp2(dtm.X,dtm.Y,dtm.z,x,y);
  
  % plotit
  
  if nargout <=  1
      h = plot3(x,y,z+5,'linewidth',1,'color','k');
  else
      pos.x = x;
      pos.y = y;
      pos.z = z;
  end