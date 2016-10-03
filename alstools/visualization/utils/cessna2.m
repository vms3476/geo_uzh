function [hp] = cessna2(x,y,z,heading,pitch,roll,scale);
% function [hp] = cessna2(x,y,z,heading,pitch,roll,scale);
% plots cessna in image using x,y,z and heading, pitch and roll and scale   
if 0
  dat = load('cessna2.dat');
  save cessna2.mat dat
else
  load cessna2.mat;
end

tris = dat(1:1117,2:3);
tris(:,1) = tris(:,1) + 1;
tris(:,2) = tris(:,2) - 1;
vert = dat(1118:1118+9680-1,:);
vert = vert*scale;
vert(:,1) = vert(:,1)+x;
vert(:,2) = vert(:,2)+y;
vert(:,3) = vert(:,3)+z;
VertNormals = dat(1118+9680:end,:);
k = 0;
for i = 1:length(tris);
  ii = [tris(i,1):tris(i,1)+tris(i,2)];
  if length(ii) <= 3
    k = k + 1;
    faces(k,:) = [ii];
  else
    for j = 1:length(ii)-2
        k = k + 1;
        faces(k,:) = ii(j:j+2);
    end
  end
end

  col = ones(size([vert(:,1),vert(:,2),vert(:,3)]))*1;
  col(:,1) = 1;
  
  
  facecolormat = ones(length(vert(:,1)),1)*[0.8 0.9 0.9];
  
  hp = patch('vertices',[vert(:,1),vert(:,2),vert(:,3)],'faces',faces,'FaceVertexCData', facecolormat, ...
             'LineStyle','none', ...
             'VertexNormals',[VertNormals(:,1),VertNormals(:,2),VertNormals(:,3)], ...
             'edgecolor','none', ...
             'facelighting','phong', ...
             'EdgeLighting','phong', ...
             'MarkerEdgeColor', 'none', ...
             'AmbientStrength', 0.5); 
  material shiny
  hold on
  xmin = find(min(vert(:,1))==vert(:,1));
  xmax = find(max(vert(:,1))==vert(:,1));
  zmax = find(max(vert(:,3))==vert(:,3));
  vertr = vert(xmin,:);
  vertg = vert(xmax,:);
  vertw = vert(zmax,:);
  ymir = find(max(vertr(:,2))==vertr(:,2));
  ymig = find(max(vertg(:,2))==vertg(:,2));
  ymiw = find(max(vertw(:,2))==vertw(:,2));
  
  % plot wingtips  
  pred = plot3(vertr(ymir,1),vertr(ymir,2),vertr(ymir,3),'.r');
  pgre = plot3(vertg(ymig,1),vertg(ymig,2),vertg(ymig,3),'.g');
  pwhi = plot3(vertw(ymiw,1),vertw(ymiw,2),vertw(ymiw,3),'.w');
  pps = [pred,pgre,pwhi];
  
  set(pps,'markersize',1,'linewidth',1);
  
  % rotate wingtips
  rotate(pps,[0 0 90],heading,[x,y,z]);
  rotate(pps,[90 0 0],pitch,[x,y,z]);
  rotate(pps,[0 90 0],roll,[x,y,z]);
  
  % rotate cessna
  rotate(hp,[0 0 90],heading,[x,y,z]);
  rotate(hp,[90 0 0],pitch,[x,y,z]);
  rotate(hp,[0 90 0],roll,[x,y,z]);

 






