% script to read in the LFI data from Rolf Meile
if 1
    cc
    lfi = xlsread('/Users/morsdorf/data/ForestGrowth/ba_peri_20140127_v2.xlsx');
    plts = unique(lfi(:,2));
    tree = unique(lfi(:,4));
    for i = 1:length(tree)
        ii = find(tree(i) == lfi(:,4));
        if length(ii) == 1
            lfitrees.plt(i) = lfi(ii,2);  
            lfitrees.x(i) = lfi(ii,5);  
            lfitrees.y(i) = lfi(ii,6);  
            lfitrees.id(i) = lfi(ii,4);  
            lfitrees.dbh(i) = lfi(ii,12);
            lfitrees.h(i) = lfi(ii,11);
            lfitrees.tim(i) = lfi(ii,3);
            lfitrees.dh(i) = NaN;
            lfitrees.ddbh(i) = NaN;
            lfitrees.dtim(i) = NaN;
            lfitrees.agr_dbh(i) = NaN;
            lfitrees.agr_h(i) = NaN;
            
        elseif length(ii) == 2
            lfitrees.plt(i) = lfi(ii(1),2);  
            lfitrees.x(i) = lfi(ii(1),5);  
            lfitrees.y(i) = lfi(ii(1),6);
            lfitrees.id(i) = lfi(ii(1),4);  
            lfitrees.h(i) = lfi(ii(1),11);
            lfitrees.dbh(i) = lfi(ii(1),12);
            lfitrees.tim(i) = complex(lfi(ii(1),3),lfi(ii(2),3));
            
            if lfi(ii(1),3) > lfi(ii(2),3) 
                lfitrees.ddbh(i) = lfi(ii(1),12) - lfi(ii(2),12);
                lfitrees.dh(i) = lfi(ii(1),11) - lfi(ii(2),11);
                lfitrees.dtim(i) = lfi(ii(1),3) - lfi(ii(2),3);
            elseif lfi(ii(1),3) < lfi(ii(2),3) 
                lfitrees.dh(i) = lfi(ii(2),11) - lfi(ii(1),11);
                lfitrees.ddbh(i) = lfi(ii(2),12) - lfi(ii(1),12);
                lfitrees.dtim(i) = lfi(ii(2),3) - lfi(ii(1),3);
            end
            lfitrees.agr_dbh(i) = lfitrees.ddbh(i)/lfitrees.dtim(i);
            lfitrees.agr_h(i) = lfitrees.dh(i)/lfitrees.dtim(i);
        elseif length(ii) > 2
            disp(['More than two inventories found for tree number ',num2str(lfi(ii(1),4)),' !'])
        end
    end
    lfitrees.h(isnan(lfitrees.h)) = 0;
    load trees2002
    alltrees02 = trees;
    
    load trees2010
    alltrees10 = trees;

    
    load dtmdsm2010.mat
    CHM = dsm;
    CHM.z = dsm.z-dtm.z;
end

k = 0;
%figure(1);clf;
%orient tall;wysiwyg;
dx = 3;
for i = 1:length(plts)
    ii = find(lfitrees.plt == plts(i));
    lfitree = subsetraw(lfitrees,ii);
    lfix = lfitree.x;
    lfiy = lfitree.y;
    cx = round(mean(lfix)/10)*10;
    cy = round(mean(lfiy)/10)*10;
    ii = alltrees02.x >= min(lfix)-dx & alltrees02.x <= max(lfix)+dx & ...
         alltrees02.y >= min(lfiy)-dx & alltrees02.y <= max(lfiy)+dx;
    jj = alltrees10.x >= min(lfix)-dx & alltrees10.x <= max(lfix)+dx & ...
         alltrees10.y >= min(lfiy)-dx & alltrees10.y <= max(lfiy)+dx;
    II = CHM.x >= cx-12-dx & CHM.x <= cx+12+dx;
    JJ = CHM.y >= cy-12-dx & CHM.y <= cy+12+dx;
    chm.x = CHM.x(II);
    chm.y = CHM.y(JJ);
    chm.z = CHM.z(JJ,II);
    if sum(ii) > 0 & sum(jj) > 0
        k = k + 1;
        %subplot(3,2,k)
        figure(k)
        subplot(1,2,1);
        myimage(chm);
        hold on
        plot(lfix,lfiy,'xb');
        ht = text(lfix,lfiy,num2str(lfitree.id'));
        set(ht,'color','r');
        %plot(alltrees02.x(ii),alltrees02.y(ii),'xg');
        %plot(alltrees10.x(jj),alltrees10.y(jj),'xb');
        axis equal;axis tight;swisstick;
        grid on
        title(['LFI Plot ',num2str(plts(i))]);
        colormap(gray);
        
        % matchtrees
        siz = 6;
        alstrees10.x = alltrees10.x(jj);
        alstrees10.y = alltrees10.y(jj);
        alstrees10.h = alltrees10.h(jj);
        alstrees02.x = alltrees02.x(ii);
        alstrees02.y = alltrees02.y(ii);
        alstrees02.h = alltrees02.h(ii);
        alstrees = matchtrees(alstrees02,alstrees10);
        alstrees.h = alstrees.alsh - alstrees.h;
        %trees02 = matchtrees(lfitree,alstrees02);
        %trees10 = matchtrees(lfitree,alstrees10);
        trees = matchtrees(lfitree,alstrees);
        plot(trees.x,trees.y,'.r','markersize',siz);
        %plot(trees10.alsx,trees10.alsy,'.m','markersize',siz);
        plot(trees.alsx,trees.alsy,'.b','markersize',siz);
        %plot(trees10.x,trees10.y,'.c','markersize',siz);
        %        plot(lfi(:,5),lfi(:,6),'xy','markersize',siz);
        %        plot(lfi(:,5),lfi(:,6),'xy','markersize',siz);
        hl = line([trees.alsx,trees.x]',[trees.alsy,trees.y]');
        set(hl,'linewidth',siz/10,'color','w');
        %hl = line([trees10.alsx,trees10.x]',[trees10.alsy,trees10.y]');
        %set(hl,'linewidth',siz/10,'color','k');
        hl = legend('LFI','ALS','location','southwest');
        subplot(1,2,2)
        %plot(trees.agr_dbh,trees.alsh,'.')
        myregress(trees.agr_dbh,trees.alsh, ...
                  {'Annual stem diameter growth rate [cm]', ...
                  'ALS estimated tree growth [m]'},[NaN,NaN]);   
        %xlabel('Annual stem diameter growth rate [cm]')
        %ylabel('ALS estimated tree growth [m]');
        grid on
        print('-depsc2',['lfi_vs_als_',num2str(plts(i)),'.eps'])
    end
end
%subplot(3,2,5)
