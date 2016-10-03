function [hcb] = pldiscret(rawf,rawl,colb)
%  function [hcb] = pldiscret(rawf,rawl,colb)
% plot first and last echo data using two colormaps (mk_powerline.m)

  if nargin == 2
      colb = 'horiz';
  end
  
  hcb = [];
  %cmap1 = ones(256,3);cmap1(:,2:3) = 0;
  %cmap2 = ones(256,3);cmap2(:,1:2) = 0;
  cmap1 = gray(256);
  cmap2 = hsv(256);
  
  dlp = [minnan(abs(rawl.z(:))),maxnan(abs(rawl.z(:)))];
  dfp = [minnan(abs(rawf.z(:))),maxnan(abs(rawf.z(:)))];
  myscatter3(rawf.x(:),rawf.y(:),rawf.z(:),rawf.z(:),cmap2,1.5);
  myscatter3(rawl.x(:),rawl.y(:),rawl.z(:),rawl.z(:),cmap1,1.5);
  colormap([cmap2;cmap1]);
  
  if ~strcmp(colb,'none')
    ax = gca;
    hcb = colorbar(colb);
    axes(hcb)
    
    llp = find_label(dlp,10);
    llp(llp < dlp(1) | llp > dlp(2)) = [];
    lfp = find_label(dfp,10);
    lfp(lfp < dfp(1) | lfp > dfp(2)) = [];
    if strcmp(colb,'vert')
      ylabel(['First Echo Hoehe [m] | Last Echo Hoehe [m]']);
      cay = get(hcb,'ylim');
      dy = (cay(2)-cay(1))/2;
      yti1 = (llp-dlp(1))/(dlp(2)-dlp(1))*dy;
      yti2 = (lfp-dfp(1))/(dfp(2)-dfp(1))*dy;
      yti2 = yti2+dy+cay(1);
      yti1 = yti1+cay(1);
      set(gca,'ytick',[yti1,yti2],'yticklabel',[llp,lfp]);
      xli = get(gca,'xlim');
      line([xli(1),xli(2)],[dy dy],'color','k','linewidth',3);
    else
      xlabel(['First Echo Hoehe [m] | Last Echo Hoehe [m]']);
      cax = get(hcb,'xlim');
      dx = (cax(2)-cax(1))/2;
      xti1 = (llp-dlp(1))/(dlp(2)-dlp(1))*dx;
      xti2 = (lfp-dfp(1))/(dfp(2)-dfp(1))*dx;
      xti2 = xti2+dx+cax(1);
      xti1 = xti1+cax(1);
      set(gca,'xaxislocation','top');
      set(gca,'xtick',[xti1,xti2],'xticklabel',[llp,lfp]);
      yli = get(gca,'ylim');
      line([yli(1),yli(2)],[dx dx],'color','k','linewidth',3);
    end
    axes(ax)
  end
  