function [] = clu2osg(fname,X,Y,Z,pl_cohull,height);
% function [] = clu2osg(fname,x,y,z,pl_cohull);
  if isempty(strfind(fname,'.osg'))
    fname = [fname,'.osg'];
  end
  if nargin >= 5
    pl_cohull = 1; % switch for computing convex hull - computes now ...
    k = 2;
  else
    pl_cohull = 0;
    k = 1;
  end
  colclass = 1; %CHANGE for single color output
  fid = fopen(fname,'W');
  
  % write header for all geometries
  
  fprintf(fid,'%s \n','Geode {');
  fprintf(fid,'%s \n','  UniqueID Geode_0 ');
  fprintf(fid,'%s \n','  DataVariance DYNAMIC ');
  fprintf(fid,'%s \n',['  name "',fname,'" ']);
  fprintf(fid,'%s \n','  cullingActive TRUE ');
  
  % only one or more trees ?
  
  if iscell(X) % more than one tree
    fprintf(fid,'%s \n',['  num_drawables ',num2str(length(X)*k),' ']);
    II = randperm(length(X));
    for i = 1:20:length(X)
      XX = []; YY = XX; ZZ = XX;
      for k = i:i+19
        mz = median(Z{k});
        ii = Z{k} > mz + 30;
        X{k}(ii)=[];Y{k}(ii)=[];Z{k}(ii)=[];
        XX = [XX;X{k}];YY = [YY;Y{k}];ZZ = [ZZ;Z{k}];
      end
      if ~isempty(XX)
        [h,x,y,z] = myscatter(XX,YY,ZZ,ZZ);
        cmap = ocean2(length(x));
        cmap2 = rand(1,3);
        for j = 1:length(x) 
          writegeometry(fid,x{j},y{j},z{j},'points',[cmap(j,:);cmap(j,:); ...
                              cmap2;cmap2]);
        end      
        if pl_cohull
          if length(X{i}) < 2000
	    prc = prctile(-Z{i},98);
	    ii = -Z{i} < prc;
	    if 1
	      [x,y,z] = convexhull(X{i}(ii)',Y{i}(ii)',Z{i}(ii)');
	      writegeometry(fid,x,y,z,'triangles',rand(1,3));  
	    else
	      [dia,vol,base,point,dum] = crdata(X{i}(ii)',Y{i}(ii)',Z{i}(ii)');
	      base = max(Z{i}(ii))-base;%+height(i);
	      [ht,hs] = treemodel(0,0,height(i),base,dia);
	      x = get(ht,'XData')';y = get(ht,'yData')';z = get(ht,'ZData')';
	      [x,y,z] = convexhull(x(:)',y(:)',z(:)');
	      z = z+min(Z{i}(ii)')-1;
	      x = x+mean(X{i}(ii)');
	      y = y+mean(Y{i}(ii)');
	      writegeometry(fid,x,y,z,'triangles',rand(1,3));  
	      x = get(hs,'XData')';y = get(hs,'yData')';z = get(hs,'ZData')';
	      [x,y,z] = convexhull(x(:)',y(:)',z(:)');
	      z = z+min(Z{i}(ii)')-1;
	      x = x+mean(X{i}(ii)');
	      y = y+mean(Y{i}(ii)');
	      writegeometry(fid,x,y,z,'triangles',rand(1,3));  
	    end
          end
        end
      end
    end
  else % one single tree but color classes
    if colclass
      X = X - minnan(X);
      Y = Y - minnan(Y);
      Z = Z - minnan(Z);
      [h,x,y,z] = myscatter(X,Y,Z,Z);
      fprintf(fid,'%s \n',['  num_drawables ',num2str(length(x))]);
      cmap = ocean2(length(x));
      for i = 1:length(x) 
	writegeometry(fid,x{i},y{i},z{i},'points',cmap(i,:));
      end      
    else
      fprintf(fid,'%s \n','  num_drawables 1 ');
      writegeometry(fid,X,Y,Z,'points');
      if pl_cohull
	[x,y,z] = convexhull(X,Y,Z);
	writegeometry(fid,x,y,z,'triangles');
      end
    end
  end
  

  fprintf(fid,'%s \n','} ');  
 
  fclose(fid)  


%-------------------------------------------------------------------------------  

function [] = writegeometry(fid,X,Y,Z,glx,col);

  if nargin == 5
    col = [1 1 1];
  end
  
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
	   '           DataVariance DYNAMIC', ...
	   '           ColorMode ON', ...
           ['           ambientColor ',num2str(col(1,:)),' 1'], ...
	   ['           diffuseColor ',num2str(col(2,:)),' 1'], ...
	   ['           specularColor ',num2str(col(3,:)),' 1'], ...
	   ['           emissionColor ',num2str(col(4,:)),' 1'], ...
	   '            shininess 0.8', ...
	   '         }',...
	   'ShadeModel {', ...
           '           DataVariance DYNAMIC', ...
           '           mode SMOOTH', ...
           '           }', ...
	   'Point      {', ...
           '           size 2', ...
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

    % Calculate Vertex Normals and smooth over faces
    
    for i = 1:3:length(X)-2
      v1 = [X(i),Y(i),Z(i)]-[X(i+1),Y(i+1),Z(i+1)];
      v2 = [X(i+2),Y(i+2),Z(i+2)]-[X(i+1),Y(i+1),Z(i+1)];
      temp = cross(v1,v2); 
      if sum(temp) == 0
	temp = [1 1 1];
      end
      temp = temp/sqrt( (temp(1)^2) + (temp(2)^2) + (temp(3)^2));
      for j = 1:3
      	C(j,i:i+2) = temp(j);
      end
    end
    [id,da] = same([X;Y;Z]','vec');
    NX = C(1,:);NY = C(2,:);NZ = C(3,:);
    for i = 1:length(id)
      if length(id{i}) > 1
	NX(id{i}) = mean(NX(id{i}));
	NY(id{i}) = mean(NY(id{i}));
	NZ(id{i}) = mean(NZ(id{i}));
      end
    end

    fprintf(fid,'%s \n',['        NormalArray  ',num2str(length(NX))]);
    fprintf(fid,'%s \n','    { ');
    for i = 1:length(NX)
      fprintf(fid,'     %6.2f %6.2f %6.2f\n',[NX(i) NY(i) NZ(i)]);
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




