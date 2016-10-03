function [] = mat3d2osg(fname,type,ax);
%function [] = mat3d2osg(fname,type,ax); 
% takes the current figure and writes the contained 'type' data to an osg file
% HINT : more than one OSG file can be loaded to osgviewer with the 
% following syntax:
% osgviewer dat1.osg dat2.osg dat3.osg ...
  
% Felix Morsdorf, RSL Zuerich, 2004-2007  

if nargin == 1
  types = {'surf','line','point','patch'};
else
  types = cellstr(type);
end

% append ending to outfilename
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
  
% loop over number of objects/colors ?

for i = 1:length(types)
  % get object information
  type = types{i};
  if nargin == 3
    hs = findobj(ax,'type',type);
  else
    hs = findobj(gcf,'type',type);
  end
  fprintf(fid,'%s \n',['  num_drawables ',num2str(length(hs)),' ']);
  
  if strcmp(lower(type),'surf') | strcmp(lower(type),'patch');
    for j = 1:length(hs) 
      cdata = get(hs(j),'cdata');
      xdata = get(hs(j),'xdata');
      ydata = get(hs(j),'ydata');
      zdata = get(hs(j),'zdata');
      if any(size(xdata) ~= size(zdata))
        [xdata,ydata] = meshgrid(xdata,ydata);
      end
      vertn = get(hs(j),'vertexnormals');
      col.face  = get(hs(j),'facecolor');
      col.facealpha  = get(hs(j),'facealpha');
      col.ambient  = col.face*get(hs(j),'ambientstrength');
      col.diffuse  = col.face*get(hs(j),'diffusestrength');
      col.specular = col.face*get(hs(j),'specularstrength');
      if strcmp(lower(type),'surf') 
        VX = squeeze(vertn(:,:,1));
        VY = squeeze(vertn(:,:,2));
        VZ = squeeze(vertn(:,:,3));  
        [m,n] = size(xdata);
        jj = [1:m-1,2:m];
        % tristripping
        xdat = reshape(xdata(jj,:),m-1,n*2);
        ydat = reshape(ydata(jj,:),m-1,n*2);
        zdat = reshape(zdata(jj,:),m-1,n*2);
        vx = reshape(VX(jj,:),m-1,n*2);
        vy = reshape(VY(jj,:),m-1,n*2);
        vz = reshape(VZ(jj,:),m-1,n*2);
      elseif strcmp(lower(type),'patch'); 
        xdat = [xdata;xdata(3)];
        ydat = [ydata;ydata(3)];
        zdat = [zdata;zdata(3)];
        vx = [vertn(1:4,1);vertn(3:4,1)];
        vy = [vertn(1:4,2);vertn(3:4,2)];
        vz = [vertn(1:4,3);vertn(3:4,3)];
        %VX = [vertn(1:4,1);vertn(1,1)];
        %VY = [vertn(1:4,2);vertn(1,2)];
        %VZ = [vertn(1:4,3);vertn(1,3)];
        xdat = xdat(:)';ydat = ydat(:)';zdat = zdat(:)';
        vx = vx(:)';vy = vy(:)';vz = vz(:)';
      end
      if min(size(xdat)) == 1
        ii = ~isnan(zdat);
        jj = ~isnan(vz);
        ii = all([ii;jj]);
        if sum(ii) > 0
          w_triangles(fid,xdat(ii),ydat(ii),zdat(ii),vx(ii),vy(ii),vz(ii),col);
        end
      else
        for I = 1:m-1
          ii = ~isnan(zdat(I,:));
          jj = ~isnan(vz(I,:));
          ii = all([ii;jj]);
          w_triangles(fid,xdat(I,ii),ydat(I,ii),zdat(I,ii),vx(I,ii),vy(I,ii),vz(I,ii),col);
        end
      end
    end
  elseif strcmp(lower(type),'line');
    for j = 1:length(hs) 
      col = get(hs(j),'color');
      x = get(hs(j),'xdata');
      y = get(hs(j),'ydata');
      z = get(hs(j),'zdata');
      if ~isempty(x)
          if strcmp(get(hs(j),'linestyle'),'none')
              siz = round(get(hs(j),'markersize')/2);
              %siz = 3;
              w_points(fid,x,y,z,col,siz);
          elseif strcmp(get(hs(j),'marker'),'none');
              siz = get(hs(j),'linewidth');
              w_lines(fid,x,y,z,col,siz);
          end
      end
    end
  end
end

fprintf(fid,'%s \n','} ');  
 
fclose(fid)  


%------------------------------------------------------------------------------

function [] = w_triangles(fid,X,Y,Z,NX,NY,NZ,col);

  if nargin == 5
    col = [1 1 1];
  end
  
  fprintf(fid,'%s \n','  Geometry { ');
  fprintf(fid,'%s \n','    DataVariance DYNAMIC ');
  fprintf(fid,'%s \n','    StateSet { ');
  fprintf(fid,'%s \n','      UniqueID StateSet_1 ');  
  fprintf(fid,'%s \n','      DataVariance STATIC ');  
  fprintf(fid,'%s \n','      rendering_hint TRANSPARENT_BIN ');
  fprintf(fid,'%s \n','      renderBinMode USE ');
  fprintf(fid,'%s \n','      GL_CULL_FACE OFF ');
  fprintf(fid,'%s \n','      GL_BLEND ON ');
  fprintf(fid,'%s \n','      GL_COLOR_MATERIAL ON ');
  fprintf(fid,'%s \n','      GL_LIGHTING OFF ');
  str = {'Material {', ...
         '           DataVariance DYNAMIC', ...
         '           ColorMode OFF', ...
         ['           ambientColor ',num2str(col.ambient),' ',num2str(col.facealpha)], ...
         ['           diffuseColor ',num2str(col.diffuse),' ',num2str(col.facealpha)], ...
         ['           specularColor ',num2str(col.specular),' ',num2str(col.facealpha)], ...
         ['           emissionColor ',num2str(col.face),' ',num2str(col.facealpha)], ...
         '            shininess 1', ...
         '         }',...
         'BlendFunc {', ...
         '   DataVariance STATIC', ...
         '   source SRC_ALPHA', ...
         '   destination ONE', ...
         ' }', ...
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
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n','    useDisplayList TRUE ');
  fprintf(fid,'%s \n','    Primitives 1 ');
  fprintf(fid,'%s \n','    { ');
  fprintf(fid,'%s \n',['    DrawArrays TRIANGLE_STRIP 0 ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n',['    VertexArray ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    { ');
  for i = 1:length(X)
    fprintf(fid,'     %6.3f %6.3f %6.3f\n',[X(i) Y(i) Z(i)]);
  end
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n','         NormalBinding PER_VERTEX ');
  fprintf(fid,'%s \n',['        NormalArray  ',num2str(length(NX))]);
  fprintf(fid,'%s \n','    { ');
  for i = 1:length(NX)
    fprintf(fid,'     %6.3f %6.3f %6.3f\n',[NX(i) NY(i) NZ(i)]);
  end
  fprintf(fid,'%s \n','    } '); 
  fprintf(fid,'%s \n','  }   ');

%------------------------------------------------------------------------------
function [] = w_points(fid,X,Y,Z,col,size);

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
  fprintf(fid,'%s \n',['      Size  ',num2str(size*5),' ']);
  fprintf(fid,'%s \n','      GL_LIGHTING ON ');
  str = {'Material {', ...
         '           DataVariance AMBIENT_DIFFUSE', ...
         '           ColorMode ON', ...
         ['           ambientColor ',num2str(col),' 1'], ...
         ['           diffuseColor ',num2str(col),' 1'], ...
         ['           specularColor ',num2str(col),' 1'], ...
         ['           emissionColor ',num2str(col),' 1'], ...
         '            shininess 1', ...
         '         }',...
         'ShadeModel {', ...
         '           DataVariance DYNAMIC', ...
         '           mode SMOOTH', ...
         '           }'};
  for i = 1:length(str)
    fprintf(fid,'%s \n',str{i});
  end
  fprintf(fid,'%s \n','      Point      {'); 
  fprintf(fid,'%s \n',['                 size ',num2str(size),' ']); 
  fprintf(fid,'%s \n','                 }'); 
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n','    useDisplayList TRUE ');
  fprintf(fid,'%s \n','    Primitives 1 ');
  fprintf(fid,'%s \n','    { ');
  fprintf(fid,'%s \n',['    DrawArrays POINTS 0 ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n',['    VertexArray ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    { ');
  for i = 1:length(X)
    fprintf(fid,'     %6.3f %6.3f %6.3f\n',[X(i) Y(i) Z(i)]);
  end
  fprintf(fid,'%s \n','    } ');         
  fprintf(fid,'%s \n','  }   ');
  
  
  
  %------------------------------------------------------------------------------
function [] = w_lines(fid,X,Y,Z,col,size);

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
  fprintf(fid,'%s \n',['      Size  ',num2str(size*5),' ']);
  fprintf(fid,'%s \n','      GL_LIGHTING ON ');
  str = {'Material {', ...
         '           DataVariance AMBIENT_DIFFUSE', ...
         '           ColorMode ON', ...
         ['           ambientColor ',num2str(col),' 1'], ...
         ['           diffuseColor ',num2str(col),' 1'], ...
         ['           specularColor ',num2str(col),' 1'], ...
         ['           emissionColor ',num2str(col),' 1'], ...
         '            shininess 1', ...
         '         }',...
         'ShadeModel {', ...
         '           DataVariance DYNAMIC', ...
         '           mode SMOOTH', ...
         '           }'};
  for i = 1:length(str)
    fprintf(fid,'%s \n',str{i});
  end
  fprintf(fid,'%s \n','      Lines      {'); 
  fprintf(fid,'%s \n',['                 width ',num2str(size),' ']); 
  fprintf(fid,'%s \n','                 }'); 
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n','    useDisplayList TRUE ');
  fprintf(fid,'%s \n','    Primitives 1 ');
  fprintf(fid,'%s \n','    { ');
  fprintf(fid,'%s \n',['    DrawArrays LINES 0 ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    } ');
  fprintf(fid,'%s \n',['    VertexArray ',num2str(length(X)),' ']);
  fprintf(fid,'%s \n','    { ');
  for i = 1:length(X)
    fprintf(fid,'     %6.3f %6.3f %6.3f\n',[X(i) Y(i) Z(i)]);
  end
  fprintf(fid,'%s \n','    } ');         
  fprintf(fid,'%s \n','  }   ');

  

%-------------------------------------------------------------------------------
function [nx,ny,nz,vx,vy,vz] = convexhull(X,Y,Z,VX,VY,VZ);
    
  dat = [X;Y;Z]';
  [k,v] = convhulln(dat);
  nx = [];ny = [];nz = [];
  vx = [];vy = [];vz = [];

  for i = 1:max(size(k))
    nx = [nx X(k(i,:))];
    ny = [ny Y(k(i,:))];
    nz = [nz Z(k(i,:))];
    vx = [vx VX(k(i,:))];
    vy = [vy VY(k(i,:))];
    vz = [vz VZ(k(i,:))];
  end 

