function [ht,XG,YG,ZG,VG] = treemodel3(x,y,z,th,ch,cd,spec,col);
%function [ht,XG,YG,ZG,VG] = treemodel3(x,y,z,th,ch,cd,spec,col);  
% draws treemodel at x,y,z, oposed to treemodel which 
% only draws in a plain of x,y ...
% treemodel3 outputs tree structure as well in voxel representation.  

% Felix Morsdorf, RSL Zuerich, Oct. 2003 - Jun. 2013
  
  % set up destination grid for voxels
  if nargout > 1  
    xyres = 0.25;
    zres = 0.25;
    sd = size(cd);
    X = min(x)-(max(cd)/2):xyres:max(x)+(max(cd)/2);
    Y = min(y)-(max(cd)/2):xyres:max(y)+(max(cd)/2);
    Z = min(z):zres:max(z)+max(th);
    [XP,YP] = meshgrid(X,Y);
    [XG,YG,ZG] = meshgrid(X,Y,Z);
    VG = zeros(size(XG));
  end
  for i = 1:length(x)
    if nargout > 1
      PG = zeros(size(XP));
    end
    [x1,y1,z1] = cylinder([0.3 0.3 0.3],4);
    x1 = x1 + x(i);
    y1 = y1 + y(i);
    z1 = z1 * (th(i)-ch(i));
    z1 = z1 + z(i);
    % tree trunk, comment out for canopy profiles !
    %[ii] = inpolygon(XP(:),YP(:),x1,y1);
    %PG(ii) = PG(ii)+1;
    %jj = find(Z >= min(z1(:)) & Z <= max(z1(:)));
    %for j = 1:length(jj)
    %  VG(:,:,jj(j)) = VG(:,:,jj(j))+PG;
    %end
    ht = 1;hf = 1;
    if nargout < 2
      ht = surf(x1,y1,z1,'Facecolor',[0.4 0.4 0],'EdgeColor','none');
      hold on
    end
    if spec(i) == 1 % paraboloid for coniferous
      env = fliplr(sqrt(0:10)/sqrt(10)*cd(i)/2);
      [x2,y2,z2] = cylinder(env,6);
      z2 = z2 * ch(i) + (th(i)-ch(i));
      z2 = z2 + z(i);
      if nargin <= 7
          facecol = [0.1 0.5 0.1];
      else
          facecol = col;
      end
    elseif spec(i) == 2 % ellipsoid for more than one crown diameter
      [x2,y2,z2] = ellipsoid(0,0,0,cd(i,1)/2,cd(i,2)/2,ch(i)/2,6);
      z2 = z2 + (th(i)-ch(i)/2);
      z2 = z2 + z(i);
      if nargin <= 7
          facecol = [0.2 0.8 0.2];
      else
          facecol = col;
      end
    end
    x2 = x2 + x(i);
    y2 = y2 + y(i);

    % voxelise tree crown information
    if nargout > 1
        for j = 1:length(Z)-1
            PG = zeros(size(XP));
            jj = find(z2(:,1) >= Z(j) & z2(:,1) <= Z(j+1));
            if ~isempty(jj)
                px = x2(jj,:);py = y2(jj,:);           
                if length(jj) > 1
                    px = mean(x2(jj,:));py = mean(y2(jj,:));
                    [ii] = inpolygon(XP(:),YP(:),px,py);
                    PG(ii) = PG(ii)+1;    
                end
            end
            VG(:,:,j) = VG(:,:,j)+PG;
        end
        ox = [x2(:);x1(:)];oy = [y2(:);y1(:)];oz = [z2(:);z1(:)];
    end
    hf = surf(x2,y2,z2,'Edgecolor','none','FaceColor',facecol,'FaceAlpha',0.5);
    if nargout == 1
      ht = [ht;hf];
      ht(i) = hf;
    end
  end
 







