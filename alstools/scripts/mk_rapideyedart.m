cc
load /Users/morsdorf/Projects/Active/3DVegLab/FabiansSims/rapideye_03sep2013.mat

load /Users/morsdorf/Projects/Active/3DVegLab/FabiansSims/dart20130903_rapideye.mat
% reik - 5,3,2
% fabian - 1,3,5

rapd(:,:,1) = flipud(rapideye_03sep2013.b5);
rapd(:,:,2) = flipud(rapideye_03sep2013.b3);
rapd(:,:,3) = flipud(rapideye_03sep2013.b2);

dart = dart20130903_rapideye(:,:,[5,3,2]);

minr = min(dart(:));
maxr = max(dart(:));
dart = (dart-minr)/(maxr-minr);
%for i = 1:3
%    dart(:,:,i) = flipud(dart(:,:,i));
%end
%minr = min(rapd(:));
%maxr = max(rapd(:));
rapd = (rapd-minr)/(maxr-minr);
for i = 1:3
    %dart(:,:,i) = adapthisteq(dart(:,:,i),'cliplimit',0.01,'nbins',1024);
    %rapd(:,:,i) = adapthisteq(rapd(:,:,i),'cliplimit',0.01,'Nbins',1024);
    %dart(:,:,i) = imadjust(dart(:,:,i));
    %rapd(:,:,i) = imadjust(rapd(:,:,i));
end
xr = rapideye_03sep2013.x;
yr = rapideye_03sep2013.y;

figure(1);clf;
subplot(1,3,1);
myimage(xr,yr,rapd);
mytitle('Rapideye @ 5m');
subplot(1,3,2);

geo = dart20130903_rapideye_info;
xd = geo.x - 2000000;
yd = geo.y - 1000000;
myimage(xd,yd,dart);
mytitle('DART @ 2m');
subplot(1,3,3);

%h = fspecial('gaussian',5,0.5);
h = fspecial('disk',3.25);
for i = 1:3
    dart(:,:,i) = imfilter(dart(:,:,i),h);
end
myimage(xd,yd,dart);
mytitle('DART @ 5m');

print -depsc2 -r600 /Users/morsdorf/Desktop/dartrapideye.eps
[XR,YR] = meshgrid(xr,yr);
[XD,YD] = meshgrid(xd,yd);
for i = 1:3
    DART(:,:,i) = interp2(XD,YD,dart(:,:,i),XR,YR);
end
return


for i = 1:3
    figure
    Y = (DART(:,:,i)*(maxr-minr))+minr;
    X = (rapd(:,:,i)*(maxr-minr))+minr;
    if 0
    % Estimate a continuous pdf from the discrete data
    [pdfx xi]= ksdensity(X(:));
    [pdfy yi]= ksdensity(Y(:));
    % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
    [xxi,yyi]     = meshgrid(xi,yi);
    [pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
    % Calculate combined pdf, under assumption of independence
    pdfxy = pdfxx.*pdfyy; 
    % Plot the results
    surf(xxi,yyi,pdfxy);
    shading interp
    set(gca,'XLim',[min(xi) max(xi)])
    set(gca,'YLim',[min(yi) max(yi)])
    view(2)
    end
    %kde2d([X(:),Y(:)]);
    myregress(double(X(:)),double(Y(:)),{'RapidEYE','DART'},[NaN NaN]);
    axis equal;
end


