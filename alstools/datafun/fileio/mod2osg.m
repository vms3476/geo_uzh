function [] = mod2osg(fname,X,Y,Z);
% function [] = mod2osg(fname,x,y,z);
% writes LIDAR model data to .osg file
% z   : length(x) by length(y) matrix
% x,y : position vectors

% Felix Morsdorf, RSL Zuerich, May 2003
  
  if isempty(strfind(fname,'.osg'))
    fname = [fname,'.osg'];
  end
  
  fid = fopen(fname,'W');
  
  % write header for all geometries
  
  fprintf(fid,'%s \n','Geode {');
  fprintf(fid,'%s \n','  UniqueID Geode_0 ');
  fprintf(fid,'%s \n','  DataVariance DYNAMIC ');
  fprintf(fid,'%s \n',['  name "',fname,'" ']);
  fprintf(fid,'%s \n','  cullingActive TRUE ');
  
  % only one or more areas
  
  if iscell(X) % more than one area
    fprintf(fid,'%s \n',['  num_drawables ',num2str(length(X)*k),' ']);
    for i = 1:length(X)
      writegeometry(fid,X{i}',Y{i}',Z{i}','triangles');
    end
  else % one single area
    fprintf(fid,'%s \n','  num_drawables 1 ');
    [m,n] = size(X);
    if min([m,n]) == 1
      [X,Y] = meshgrid(X,Y);
      [m,n] = size(X);
    end
    
    NX=[];NY=[];NZ=[];
    for j = 1:n-1
      for i = 1:m-1
      	NX =  [NX X(i,j) X(i+1,j) X(i,j+1) X(i,j+1) X(i+1,j) X(i+1,j+1)];
	NY =  [NY Y(i,j) Y(i+1,j) Y(i,j+1) Y(i,j+1) Y(i+1,j) Y(i+1,j+1)];
	NZ =  [NZ Z(i,j) Z(i+1,j) Z(i,j+1) Z(i,j+1) Z(i+1,j) Z(i+1,j+1)];
      end 
    end
    writegeometry(fid,NX(:)',NY(:)',NZ(:)','triangles');
  end

  fprintf(fid,'%s \n','} ');  
 
  fclose(fid)  


%-------------------------------------------------------------------------------  

function [] = writegeometry(fid,X,Y,Z,glx);

  fprintf(fid,'%s \n','  Geometry { ');
  fprintf(fid,'%s \n','    DataVariance DYNAMIC ');
  fprintf(fid,'%s \n','    StateSet { ');
  fprintf(fid,'%s \n','      UniqueID StateSet_1 ');  
  fprintf(fid,'%s \n','      DataVariance STATIC ');  
  fprintf(fid,'%s \n','      rendering_hint OPAQUE_BIN ');
  fprintf(fid,'%s \n','      renderBinMode INHERIT ');
  fprintf(fid,'%s \n','      GL_CULL_FACE OFF ');
  fprintf(fid,'%s \n','      GL_LIGHTING ON ');
  if 1%strcmp(glx,'triangles')
    str = {'Material {', ...
	   '           DataVariance STATIC', ...
	   '           ColorMode ON', ...
	   ['           ambientColor ',num2str(rand(1,3)),' 1'], ...
	   ['           diffuseColor ',num2str(rand(1,3)),' 1'], ...
	   ['           specularColor ',num2str(rand(1,3)/2),' 1'], ...
	   ['           emissionColor ',num2str(rand(1,3)/2),' 1'], ...
	   '            shininess 0', ...
	   '         }',...
	   'ShadeModel {', ...
           '           DataVariance STATIC', ...
           '           mode SMOOTH', ...
           '           }', ...
	   'Point      {', ...
           '           size 1', ...
           '           }'};
    for i = 1:length(str)
      fprintf(fid,'%s \n',str{i});
    end
  end
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n','    useDisplayList TRUE ');
  fprintf(fid,'%s \n','    Primitives 1 ');
  fprintf(fid,'%s \n','    { ');

  switch glx
    case 'points'
     fprintf(fid,'%s \n',['    DrawArrays POINTS 0 ',num2str(length(X)),' ']);
    case 'triangles'
     fprintf(fid,'%s \n',['    DrawArrays TRIANGLES 0 ',num2str(length(X)),' ']);
  end
  
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n',['    VertexArray ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    { ');
  for i = 1:length(X)
    fprintf(fid,'     %8.3f %8.3f %8.3f\n',[X(i) Y(i) Z(i)]);
  end
  fprintf(fid,'%s \n','    } ');
  if strcmp(glx,'triangles')
    fprintf(fid,'%s \n','         NormalBinding PER_VERTEX ');
    fprintf(fid,'%s \n',['        NormalArray  ',num2str(length(X))]);
    fprintf(fid,'%s \n','    { ');

    % Calculate Vertex Normals and smooth over faces
    
    for i = 1:3:length(X)-2
      v1 = [X(i),Y(i),Z(i)]-[X(i+1),Y(i+1),Z(i+1)];
      v2 = [X(i+2),Y(i+2),Z(i+2)]-[X(i+1),Y(i+1),Z(i+1)];
      temp = cross(v1,v2);     
      temp = temp/sqrt( (temp(1)^2) + (temp(2)^2) + (temp(3)^2));
      for j = 1:3
      	C(j,i:i+2) = temp(j);
      end
    end
    [id,da] = same([X;Y;Z]','vec');
    NX = C(1,:);NY = C(2,:);NZ = C(3,:);
    for i = 1:length(id)
      NX(id{i}) = -mean(NX(id{i}));
      NY(id{i}) = -mean(NY(id{i}));
      NZ(id{i}) = -mean(NZ(id{i}));
    end

    for i = 1:length(X)-2
      fprintf(fid,'     %8.3f %8.3f %8.3f\n',[NX(i) NY(i) NZ(i)]);
    end
    fprintf(fid,'%s \n','    } ');
  end
  
  fprintf(fid,'%s \n','  }   ');


%-------------------------------------------------------------------------------
function [nx,ny,nz] = convexhull(X,Y,Z);
    
  dat = [X;Y;Z]';
  [k,v] = convhulln(dat);
  nx = [];ny = [];nz = [];

  for i = 1:max(size(k))
    nx = [nx X(k(i,:))];
    ny = [ny Y(k(i,:))];
    nz = [nz Z(k(i,:))];
  end 




