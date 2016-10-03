function [hg] = swisstick(lang);
% Makes nice x/y ticks and labels for swiss national 
% coordinate system  
  
  xti = get(gca,'xtick');
  yti = get(gca,'ytick');
  set(gca,'xticklabel',[],'yticklabel',[]);
  sx = fix(min(xti)/1000)*1000;
  sy = fix(min(yti)/1000)*1000;
  xti = xti - sx;
  yti = yti - sy;
  set(gca,'xticklabel',xti);  
  set(gca,'yticklabel',yti);  
  if nargin == 0 
      xlabel(['Easting (-',num2str(sx),' m)'])
      ylabel(['Northing (-',num2str(sy),' m)'])
  else
      xlabel(['Rechtswert (-',num2str(sx),' m)'])
      ylabel(['Hoehe (-',num2str(sy),' m)'])
  end
  
  