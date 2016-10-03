% script for animating the sun azimuth and elevation over a dtm
tt.year = 2007;
tt.month = 6;
tt.day = 21;
tt.hour = 7;
tt.min = 1;
tt.sec = 0;
tt.UTC = +1;
%Ofenpass gross
loc.latitude = 46.6592;
loc.longitude = 10.2346;
loc.altitude = 1887;

%Hinwil Small
%loc.latitude = 47.297140;
%loc.longitude = 8.824373;
%loc.altitude = 555;
hour = 4:22;
minu = 0:5:55;
[hour,minu] = meshgrid(hour,minu);
hour = hour(:);
minu = minu(:);
pos = get(gca,'position');
axp = gca;
hat = axes('position',pos,'visible','off','units','normalized');
for i = 1:length(hour)
  tt.hour = hour(i);
  tt.min = minu(i);
  sun = sun_position(tt,loc);
  az(i) = sun.azimuth;
  el(i) = sun.zenith;
  hl = lightangle(180-sun.azimuth,90-sun.zenith);
  axes(hat)
  ht(1) = text(0.05,0.25,['Azimuth   : ',num2str(az(i))]);
  ht(2) = text(0.05,0.2,['Elevation : ',num2str(90-el(i))]);
  ht(3) = text(0.05,0.15,['Time      : ',num2str(hour(i)),':', ...
                      int2strv(minu(i),1,'0',2)]);
  set(ht,'fontweight','bold','fontsize',12);
  drawnow
  eval(['print -djpeg95 /Users/morsdorf/Desktop/sunani/ani_', ...
        int2strv(i,1,'0',4),'.jpg'])
  delete(hl);
  delete(ht)
  axes(axp)
end
