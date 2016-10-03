function [hp,hcb] = mycatter2(X,Y,C,cmap,siz);
% [h] = fscatter2(X,Y,C,cmap,siz);
% Plots point cloud data in cmap color classes and 3 Dimensions,
% much faster and very little memory usage compared to scatter3 !
% X,Y,Z,C are vectors of the same length
% with C being used as index into colormap (can be any values though)
% cmap is optional colourmap to be used
% h are handles to the line objects

% Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich

% Lukas, Sept 04, some minor changes, makes it a bit faster (->profiler)

  if nargin <= 2
    if nargin == 1
      cmap = ocean2(256);
      numclass = 256;
    else
      cmap = Y;
      numclass = length(Y);
    end
    if isfield(X,'int')
      C = X.int;
      Z = X.z;
    else
      C = X.z;
      Z = C;
    end
    Y = X.y;
    NX = X.x;
    clear X;
    X = NX;clear NX;
    siz = 5;
  elseif nargin == 4
    numclass = 256; % Number of color classes
    cmap = hsv(256);
    siz = 5;
  elseif nargin == 5
    numclass = max(size(cmap));
    siz = 5;
    if numclass == 1
      siz = cmap;
      cmap = hsv(256);
      numclass = 256;
    end  
  elseif nargin == 6
    numclass = max(size(cmap));
    if numclass == 1
      siz = cmap;
      cmap = hsv(256);
      numclass = 256;
    end  
  end
  
  % avoid too many calculations
  
  mins = min(C);
  maxs = max(C);
  minz = min(Z);
  maxz = max(Z);
  minx = min(X);
  maxx = max(X);
  miny = min(Y);
  maxy = max(Y);
  
  % construct colormap :
  
  col = cmap;
  
  % determine index into colormap (viel schneller, müsste aber genauso
  % effektiv sein? Klasse in nummclass aufgeteilt. Weiss nicht so genau was
  % interp macht)
  
  ii = floor( (C - mins ) * (numclass-1) / (maxs - mins) );
  ii = ii + 1;
  %old code--------------------------------------------------------
  %ii = round(interp1([floor(mins) ceil(maxs)],[1 numclass],C));
  %------------------------------------------------------------------
  hold on
  colormap(cmap);
  %axis([minx maxx miny maxy minz maxz]);

  
  % plot each color class in a loop
  %old code-----------------------------------------------------
  %k = 0;
  %for j = 1:numclass
  %  jj = find(ii == j);
  %  if ~isempty(jj)
  %    k = k + 1;
  %    h(k) = plot3(X(jj),Y(jj),Z(jj),'.','color',col(j,:), ...
  %		 'markersize',.5);
  %end  
  %end
  %--------------------------------------------------------------
  
  %plot dauert etwas länger, jedoch wegen fehlender find-funktion insgesammt
  %kürzer (->profiler)
  
  hold on
  k = 0;o = k;
  for j = 1:numclass
    jj = (ii(:)== j);
    if ~isempty(jj)
      k = k + 1;
      h = plot3(X(jj),Y(jj),Z(jj),'.','color',col(j,:),'markersize',siz);
      if ~isempty(h)
        o = o+1;
        hp(o) = h;
      end
    end  
  end
  caxis([min(C) max(C)]);
  axis image;rotate3d on;view(3);
  box on;
  %  hcb = colorbar('location','east');


