function [] = animod(fname,area);
% function [] = animod(fname,area);
%   animate laser model with dsm as texture

if nargin == 0
   load /data/sar53/morsdorf/ilfuorn3d.mat
   load /data/sar53/morsdorf/ilfuornrgb.mat
   rgb = double(ilfuornrgb);
   % scale rgb
   for i = 1:3
     rgb(:,:,i) = rgb(:,:,i)/max(max(rgb(:,:,i)));
     rgb(:,:,i) = histeq(rgb(:,:,i));
   end
   area = [ min(X) max(X) min(Y) max(Y)];
end

[X,Y,Z]=ncread(fname, ...
		 {'xdist','ydist','z'}, ...
		 'xdist',area(1:2),'ydist',area(3:4));
fname = strrep(fname,'model/dtm','rgb');
fname = strrep(fname,'DTM','RGB');

[rgb,x,y] = rgbload(fname,area);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5; 
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
dz = 10;
h = fspecial('average',dz);
MZ = imfilter(Z,h,'symmetric');
%DZ = Z;%-MZ;
[IX,IY] = meshgrid(x,y);
%keyboard
[NX,NY] = meshgrid(linspace(min(X),max(X),length(X)*2), ...
		   linspace(max(Y),min(Y),length(Y)*2));
%keyboard
if 0
[C,S] = wavedec2(Z,5,'db3');
co = 50;
pc = sign(C);
%NC(abs(C) < 40) = 0;  
[NC,idx] = sort(abs(C));
[dum,idx2] = sort(idx); %Resort coefficients
NC(1:end-1000) = 0;
NC = NC(idx2).*pc;
%[NC]  = wthcoef2('h',C, S,[2 3 4 5],[co co co co],'s');
%[NC]  = wthcoef2('v',NC,S,[2 3 4 5],[co co co co],'s');
%[NC]  = wthcoef2('d',NC,S,[2 3 4 5],[co co co co],'s');
%[WZ]  = waverec2(NC,S,'db3');

end

for i = 1:3
  DZ(:,:,i) = interp2(IX,IY,squeeze(rgb(:,:,i)),NX,NY);
end


[DZ,map] = rgb2ind(DZ,256,'nodither');
props.Cdata = (DZ);
%X = flipud(X);
surface(X(1:dz:end),Y(1:dz:end),MZ(1:dz:end,1:dz:end)',props);
light('position',[-1 0 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);
view(3)

grid on
set(gcf,'renderer','opengl')
set(gca,'visible','off')
set(gcf,'units','normalized');
set(gca,'position',[0 0 1 1]);
material dull
colormap(map)
%caxis([min(min(DZ)) max(max(DZ))])
brighten(-0.3)
%view(0,90)
%for i = 1:500;camorbit(1,0);drawnow;end
