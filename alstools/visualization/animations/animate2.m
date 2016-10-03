% animates 3d window throuh rotation and saves data as single tiffs
  ix = 2048;iy = 1024;
  %set(gcf,'position',[0.1 0.1 0.8 0.8]);
  set(gcf,'position',[1 1 ix iy])
  set(gcf,'paperunits','centimeters','paperposition',[1 1 ix/50 iy/50],'papersize',[ix/50 iy/50]);
  %orient landscape
  wysiwyg
  %axis off;axis tight
  %orient landscape;
  set(gcf,'color',[0.5 0.5 0.5],'inverthardcopy','off')
  for J = 1:2
    subplot(1,2,J)
    axis vis3d;
    xlabel('');ylabel('');zlabel('');
  end
  %wysiwyg
  %h = waitbar(0,'Please grab a cup of coffee ... ');
  az = linspace(-360,360,1800);
  el = [ones(1,300)*30,linspace(30,0,300),linspace(0,70,600),linspace(70,90,300),linspace(90,30,300)];
  alph = sawtooth(linspace(0,8*pi,1800),0.5)+1;
  alph( alph > 1 ) = 1;
  alph = 1-alph;
  for I = 1:length(az);
    for J = 1:2
      subplot(1,2,J)
      view(az(I),el(I));
      %set(ht,'facealpha',alph(I),'edgealpha',alph(I));
    end
    eval(['print -dtiff ani_mod',int2strv(I,1,'0',4),'.tif'])
  end
  %close(h)
