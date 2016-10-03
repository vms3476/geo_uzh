cd /Users/morsdorf/Desktop/TLS_Laegern/
filez = dir('*.las');
for i = 1:length(filez);
    [s,h,v] = readlas(filez(i).name,'xyzRGB');
    clf;
    subs = 1;
    fscatterRGB(s.X(1:subs:end),s.Y(1:subs:end),s.Z(1:subs:end),s.r(1:subs:end),s.g(1:subs:end),s.b(1:subs:end));
    drawnow;
    mat3d2osg(filez(i).name)
    if i == 1 
            subs = 4;clf;
            fscatterRGB(s.X(1:subs:end),s.Y(1:subs:end),s.Z(1:subs:end),s.r(1:subs:end),s.g(1:subs:end),s.b(1:subs:end));
            drawnow;
            mat3d2osg([filez(i).name,'thinned'])
    end
end
