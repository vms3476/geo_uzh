function [h] = shadmod(x,y,z,c);
% function [h] = shadmod(x,y,z,c);
% plots data model in shaded way
  if nargin == 3
    c = z;
  elseif nargin == 1
    y = x.y;
    z = x.z;
    x = x.x;
    c = z;
  end
  

props.AmbientStrength = 0.5;
props.DiffuseStrength = 0.4;
props.SpecularColorReflectance = 1; 
props.SpecularExponent = 2;
props.SpecularStrength = 0.7;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.cdata = c;
h = surface(x,y,z,props);

view(0,90);
%colormap(myspecmap)
%colormap(gray)
if min(size(c)) == 1
   colormap(ocean2)
end
set(gca,'color','k','xlim',[min(x(:)) max(x(:))], ...
        'ylim',[min(y(:)) max(y(:))])

%swisstick;
dasp = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[1 1 0.11]);
axis tight
material dull;
hl = light;
lightangle(hl,225,60);
h = [h,hl];

