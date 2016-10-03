% script for producing the reference coordinates of berner jura
cc
cd /Volumes/Data3/BernerJura
str = {'Bevilard','Cortebert','Moutier','Orvin','Vauffelin'};
for i = 1:length(str);
    figure(i);clf
    if 0
        % load point data 
        fpath = '/Volumes/Data3/BernerJura/Referenz/Kontrollflaechen/';
        fid = fopen([fpath,str{i},'/',str{i},'.csv'],'r');
        dat = textscan(fid,'%s%f%f%f%f%s', 'Delimiter',';');
        fclose(fid);
        % get area, including buffer
        buf = 5;
        are = [min(dat{3})-buf max(dat{3})+buf min(dat{4})-buf max(dat{4})+buf];
        file = [num2str(floor((are(1)-2000000)/1000)),'_',num2str(floor((are(3)-1000000)/1000)), '.las'];
        disp(['Loading raw data from file :',file]);
        raw = readlas(['/Volumes/Data3/BernerJura/Punktwolke_klassiert/',file]);
        ii = raw.x >= are(1) & raw.x <= are(2) & raw.y >= are(3) & raw.y <= are(4);
        raw = subsetraw(raw,ii);
        cd /Volumes/Data3/BernerJura/DOM/ArcInfo
        [bbox,fname] = bboxasciigrid(dir('*.asc'),are);
        disp(['Loading DOM data from file :',fname]);
        dom = asciigrid2mat(fname);
        cd /Volumes/Data3/BernerJura/DTM/ArcInfo
        [bbox,fname] = bboxasciigrid(dir('*.asc'),are);
        disp(['Loading DTM data from file :',fname]);
        dtm = asciigrid2mat(fname);
        ii = find(dom.x >= are(1) & dom.x <= are(2));
        jj = find(dom.y >= are(3) & dom.y <= are(4));
        DOM.z = dom.z(jj,ii);
        DOM.x = dom.x(ii);
        DOM.y = dom.y(jj);
        [DOM.X,DOM.Y] = meshgrid(DOM.x,DOM.y);
        DTM.z = dtm.z(jj,ii);
        DTM.x = dtm.x(ii);
        DTM.y = dtm.y(jj);
        [DTM.X,DTM.Y] = meshgrid(DTM.x,DTM.y);
        clear dom dtm
        eval(['save /Volumes/Data3/BernerJura/refpoints_',str{i},'.mat DOM DTM raw dat ' ...
              ' are']);
    else
        load(['/Volumes/Data3/BernerJura/refpoints_',str{i},'.mat']);   
    end
    shadmod(DOM);
    hold on
    pz = interp2(DOM.X,DOM.Y,DOM.z,dat{3},dat{4});
    MEX(i) = mean(dat{3});
    MEY(i) = mean(dat{4});
    plot3(dat{3},dat{4},dat{5},'ow');
    plot3(dat{3},dat{4},dat{5},'xk');
    text(dat{3},dat{4},dat{5}+1,[num2str(dat{2})]);
    title(['Hausdach in Umgebung von ',str{i}]);
    swisstick
    eval(['print -depsc2 -r600 DSM_Haus_Referenz_',str{i},'_topview.eps'])
    axis equal;axis tight;
    if 1
        view(14,40);
        eval(['print -depsc2 -r600 DSM_Haus_Referenz_',str{i},'_sideview.eps'])
        clf
        
        if i == 1
            planes = {[2021,2022,2020],[2023,2022,2020],[2020,2003,2023],[2023, ...
                                2003,2002]};
        elseif i == 2
            planes = {[2020,2006,2004],[2020,2018,2004],[2022,2020,2008],[2020, ...
                            2008,2006]};
        elseif i == 3
            planes = {[2014,2002,2032],[2002,2003,2032],[2033,2014,2001],[2014, ...
                                2002,2001]};
        elseif i == 4
            planes = {[2010,2006,2004],[2006,2004,2005],[2009,2008,2002],[2008, ...
                                2007,2002]};
        elseif i == 5
            planes = {[2014,2002,2013],[2013,2002,2003],[2014,2002,2001],[2014, ...
                                2010,2001]};
        end
        viewangs = [18.5 0;-5.5 0; -0.5 0; -78.5 0;-58.5 0];
        mex = mean(raw.x);
        mey = mean(raw.y);
        orient landscape;wysiwyg;
        subplot(1,3,1)
        myscatter3(raw.x-mex,raw.y-mey,raw.z,raw.z);
        hold on
        plot3(dat{3}-mex,dat{4}-mey,dat{5},'ok','markersize',3);
        cols = jet(4);
        for k = 1:length(planes)
            for l = 1:3
                ii = find(dat{2} == planes{k}(l));
                X(l) =  dat{3}(ii);
                Y(l) =  dat{4}(ii);
                Z(l) =  dat{5}(ii);
            end
            hp = patch(X-mex,Y-mey,Z,cols(k,:),'facealpha',0.5);
            PLANES{k} = [X',Y',Z'];
        end
        axis equal;axis tight;
        view(2)
        hcb = colorbar('South');
        pos = get(hcb,'position');
        
        %eval(['print -depsc2 -r600 RAW_Haus_Referenz_',str{i},'_topview.eps'])
        ax1 = gca;
        subplot(1,3,2)
        ax2 = gca;
        c = copyobj(allchild(ax1),ax2);
        axis equal;axis tight;
        view(viewangs(i,1),viewangs(i,2));
        title(str{i})
        hcb2 = colorbar('SouthOutside');
        pos = get(hcb2,'position');
        delete(hcb2);
        subplot(1,3,3);
        ax3 = gca;
        c = copyobj(allchild(ax1),ax3);
        axis equal;axis tight;
        view(3)
        %    set(hcb,'position',[pos(1)*3 pos(2)*0.94 pos(3) pos(4)*0.5]);
        pos(2) = pos(2)*2.7;
        pos(4) = pos(4)*0.5;
        set(hcb,'position',pos);
        drawnow;
        %pos = get(hcb,'position');
        %set(hcb,'position',[pos(1)*1.5 pos(2)*0.7 pos(3) pos(4)]);
        %eval(['print -depsc2 -r600 RAW_Haus_Referenz_',str{i},'_sideview.eps'])
        eval(['print -depsc2 -r600 RAW_Haus_Referenz_',str{i},'_threeviews.eps'])
        
        %mat3d2osg(['Haus_',str{i}]);
        %clf
        %raw.Classes = raw.Classification;
        %raw.Classes(raw.Classification == 2) = 1;
        %raw.Classes((raw.Classification >= 3 & raw.Classification <= 5)) = 2;
        %raw.Classes(raw.Classification == 6) = 3;
        %raw.Classes(raw.Classification == 0 | raw.Classification > 6) = 4;
        %col = [0.7 0.7 0;0 1 0;1 0 0;0.3 0.3 0.3]; 
        %for k = 1:4
        %    ii = raw.Classes == k; 
        %    hp(k) = plot3(raw.x(ii)-mex,raw.y(ii)-mey,raw.z(ii),'.');
        %    set(hp(k),'color',col(k,:),'markersize',10);
        %    hold on
        %end
        %axis equal;axis tight
        %legend(hp,'Boden','Vegetation','Gebaeude','Andere','location','northeast');
        %mat3d2osg(['Haus_classes_',str{i}]);
        %title(['Klassifikationsergebnis bei ',str{i}]);
        %eval(['print -depsc2 -r600 RAW_Haus_Referenz_',str{i}, ...
        %      '_classification.eps']);


    stats = coregpl(raw,PLANES);
    eval([str{i},'.stats = stats;']);
    %ii = [];
    %for k = 1:length(stats)
    %    if stats(k).R(1,1) > 0.99
    %        ii = [ii,k];
    %    end
    %end
    stpp(i)  = mean([stats(:).std]);
    meanp(i) = mean([stats(:).mean]);
    diffp = mean([stats(:).t]');
    if i == 2
        diffp = mean([stats([1,2,4]).t]');
    end
    results(i,:) = [diffp]*100;
    end
end
fid = fopen('/Volumes/Data3/BernerJura/results_match.csv','w');
for i = 1:length(results)
    fprintf(fid,'%2.6f,%2.6f,%2.6f\n',[results(i,:)]);
end
fclose(fid);




