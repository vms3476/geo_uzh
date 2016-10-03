close all
[optech]=imread('/Users/morsdorf/Desktop/optech_zermatt.tif');
optech = flipud(squeeze(optech(:,:,1)));
[riegl]=imread('/Users/morsdorf/Desktop/riegl_zermatt.tif');
riegl = flipud(squeeze(riegl(1:end-1,1:end-1,1)));

NDI = (riegl)./(optech);


figure;orient tall;wysiwyg;
[m,n] = size(NDI);
x = linspace(624000, 624500, n);
y = linspace(96500,97000, m);
subplot(3,1,1)
myimage(x,y,optech);
colormap(gray);
swisstick
grid on
title('Optech ALTM Gemini @ 1064 nm');
hcb = colorbar;
axes(hcb);
ylabel('Intensity []');

subplot(3,1,2)
myimage(x,y,riegl);
colormap(gray);
swisstick; grid on
%xlabel('');
title('Riegl LMS Q560 @ 1550 nm');
hcb = colorbar;
axes(hcb);
ylabel('Intensity []');


subplot(3,1,3)
myimage(x,y,NDI);
caxis([0 5])
colormap(gray);
swisstick; grid on
%xlabel('');
title('Simple ratio 1550/1064');
hcb = colorbar;
axes(hcb);
ylabel('Intensity ratio []');
print -dpdf -r600 /Users/morsdorf/Publications/IEEE_MSLZermatt/mint_sr_zermatt_bw.pdf
colormap(hsv)
print -dpdf -r600 /Users/morsdorf/Publications/IEEE_MSLZermatt/mint_sr_zermatt_col.pdf
