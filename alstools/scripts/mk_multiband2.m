% script for generating multistack images for ALS data in cdf format

res = 2;numlay = [0,0.5,1,1.5,2,4,8,16,32];  
dtm = loadmodel('/data/hinwil_big/model/dtm/Hinwil_DTM.cdf'); 
if ~isfield(dtm,'X')
      [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
end

xg = floor(min(dtm.x)):res:ceil(max(dtm.x));
yg = floor(min(dtm.y)):res:ceil(max(dtm.y));
[XG,YG] = meshgrid(xg,yg);
[m,n] = size(XG);

% boxes for processing 
xb = [min(xg):1000:max(xg),max(xg)];
yb = [min(yg):1000:max(yg),max(yg)];

if length(numlay)>1
  numl = length(numlay);
else
  numl = numlay;
end

omultih = ones(m,n,numl)*NaN;
omultid = ones(m,n,numl)*NaN;
onumhits = ones(m,n)*NaN;

for I = 1:length(xb)-1
  for J = 1:length(yb)-1
    [raw.x,raw.y,raw.z] = getrawdata('/data/hinwil_big/raw',[xb(I) xb(I+1) yb(J) yb(J+1)]);
    raw.z = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
    ii = xg >= min(raw.x) & xg <= max(raw.x);
    jj = yg >= min(raw.y) & yg <= max(raw.y);
    sdtm.x = xg(ii);
    sdtm.y = yg(jj);
    [x,y,multih,multid,numhits] = als2multi(sdtm,raw,-1,numlay);
    hw = waitbar(0,'Sorting elements to raster ...');
    [mm,nn] = size(multih(:,:,1));
    for i = 1:nn
      ii = x(i) == xg;              
      for j = 1:mm
        jj = y(j) == yg;
        if sum(ii) == 1 & sum(jj) == 1
          if ~isnan(multih(j,i,:))
            omultih(jj,ii,:) = multih(j,i,:);
          end
          if ~isnan(multid(j,i,:))
            omultid(jj,ii,:) = multid(j,i,:);
          end
          if ~isnan(numhits(j,i))
            onumhits(jj,ii) = numhits(j,i);
          end
        end
        if fix(i*j/1000) == (i*j/1000)
          waitbar((i*j)/(mm*nn),hw);
        end
      end
    end
    close(hw);
  end
end
save multimosaic_hinwil_toposys.mat xg yg omultih omultid onumhits 


