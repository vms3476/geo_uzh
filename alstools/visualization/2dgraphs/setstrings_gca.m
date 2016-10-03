function [] = setstrings(fs,col);
%function [] = setstrings(fs,col);
  if nargin == 0
    fs = 12;
    col = 'k';
  elseif nargin == 1
    col = 'k';
  end
    
  fontname = 'Helvetica';
  hstr = findobj(gca,'type','text');
  set(hstr,'fontweight','bold','fontsize',fs,'fontname',fontname,'color',col);
  hstr = findobj(gca,'type','string');
  set(hstr,'fontweight','bold','fontsize',fs,'fontname',fontname,'color',col);
  hgca = findobj(gca,'type','axes');
  set(hgca,'fontweight','bold','fontsize',fs,'fontname',fontname);
  hgcay = get(hgca,'ylabel');
  hgcax = get(hgca,'xlabel');
  hgcaz = get(hgca,'zlabel');
  hg = [];
  if iscell(hgcay)
    for i = 1:length(hgcay)
      hg = [hg hgcay{i}];
    end
  else
    hg = hgcay;
  end
  if iscell(hgcax)
    for i = 1:length(hgcax)
      hg = [hg hgcax{i}];
    end   
  else
    hg = [hg,hgcax];
  end
  if iscell(hgcaz)
    for i = 1:length(hgcaz)
      hg = [hg hgcaz{i}];
    end   
  else
    hg = [hg,hgcaz];
  end
  set(hg,'fontweight','bold','fontsize',fs,'fontname',fontname);





