function [hp,hcb] = fscatterRGB(X,Y,Z,r,g,b,siz);
% [hp,hcb] = fscatterRGB(X,Y,Z,R,G,B,siz);
% Plots point cloud data in cmap color classes and 3 Dimensions,
% much faster and very little memory usage compared to scatter3 !
% X,Y,Z,C are vectors of the same length
% with C being used as index into colormap (can be any values though)
% cmap is optional colourmap to be used
% h are handles to the line objects

% Felix Morsdorf, Jan 2003 - Apr 2015, Remote Sensing Laboratory, Zurich

% Lukas, Sept 04, some minor changes, makes it a bit faster (->profiler)
  numclass = 4096;
  if nargin <= 2
      Y = X.Y;
      Z = Z.Z;
      r = X.r;
      g = X.g;
      b = X.b;
      X = X.X;
      clear X;
      siz = 3;
  elseif nargin == 6
      siz = 3;
  end
  [IND,cmap] = rgb2ind(reshape([r g b],length(r),1,3),numclass);

  % avoid too many calculations
  
  % construct colormap :
  
  col = cmap;
  
  % determine index into colormap (viel schneller, müsste aber genauso
  % effektiv sein? Klasse in nummclass aufgeteilt. Weiss nicht so genau was
  % interp macht)

  hold on
  colormap(cmap);
  IND = IND + 1;
  hold on
  k = 0;
  for j = 1:numclass
    jj = (IND(:)== j);
    if ~isempty(jj)
      h = plot3(X(jj),Y(jj),Z(jj),'.','color',col(j,:),'markersize',siz);
      if ~isempty(h)
        k = k + 1;
        hp(k) = h;
      end
    end  
  end

  axis image;rotate3d on;view(3);
  box on;



