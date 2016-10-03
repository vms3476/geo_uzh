clf
cd /Users/morsdorf/Desktop/APEX_PCA
img = enviread('pcaOut7b_20140718_geo');
A = squeeze(img.z(:,:,[1,3,4]));

for i = 1:3
    mina = min(min(A(:,:,i)));
    maxa = max(max(A(:,:,i)));
    A(:,:,i) = (A(:,:,i)-mina)/(maxa-mina);
    A(:,:,i) = histeq(A(:,:,i),256);
end
are = [669500 674500 258500 260250];
orient landscape
%wysiwyg;
subplot(2,1,2);
myimage(img.x-2000000,img.y-1000000,A);
axis(are);
swisstick;
mytitle('IS - APEX PCA : 1 <> R, 3 <> G, 4 <> B');
%print -depsc2 -r600 apex_pca_134.eps

%mhe = enviread('/Users/morsdorf/data/Laegeren_Multiband_ALS/multi_height_laegeren');

mra = enviread('/Users/morsdorf/data/Laegeren_Multiband_ALS/multi_ras_laegeren');


A = squeeze(mra.z(:,:,[1,2,3]));

for i = 1:3
    mina = min(min(A(:,:,i)));
    maxa = max(max(A(:,:,i)));
    A(:,:,i) = (A(:,:,i)-mina)/(maxa-mina);
    A(:,:,i) = histeq(A(:,:,i),256);
end

subplot(2,1,1);
myimage(mra.x,mra.y,A);
axis(are);
swisstick;
mytitle('ALS : mean vegetation height <> R, std. vegetation height <> G, mean intensity <> B');
%set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))) 
print -dpdf -r900 apex_pca_134_als_MHSHMI.pdf