if 1 
    cd /Volumes/lidarlab/data/KantonAargau/Mat_Files_comb/
    filezd = dir('*.mat');
    cd /Volumes/lidarlab/data/KantonAargau/Mat_Files_Statistics/raster_statistics
    filez = dir('*.mat');
    RGB = ones(1000,1000,3);
    info = geotiffinfo('DTM_622000_264000.tif');
    for i = 432:length(filez)
        %nm = filez(i).name;
        %east = str2num(nm(1:6));
        %nort = str2num(nm(8:13));
        %if east == 635000 & nort == 264000
            load(filez(i).name);
            X = structure.x_ras_S;
            Y = structure.y_ras_S;
            RGB(:,:,1) = (0.8*structure.veg_height_mean_S./40)+0.2;
            RGB(:,:,2) = (0.8*structure.veg_height_std_S./30)+0.2;
            RGB(:,:,3) = structure.int_mean_S/200;
            RGB(RGB>1)=1;
            RGB(RGB<0)=0;
            RGB = RGB*255;
            RGB(isnan(RGB)) = 0;
            %for j = 1:3;
            %    RGB(:,:,j) = histeq(RGB(:,:,j));
            %end
            
            %figure(1);clf;
            %myimage(X,Y,RGB);
            %drawnow;
            if 0
                load(['/Volumes/lidarlab/data/KantonAargau/Mat_Files_comb/',filezd(i).name]);
                DTM.x=data.dtm_x(1:2:end);
                DTM.y=data.dtm_y(1:2:end);
                DTM.Z=imresize(data.dtm,[1000,1000]);
                
                figure(1);clf;
                [DTM.X,DTM.Y] = meshgrid(DTM.x,DTM.y);
                shadmod(DTM.X,DTM.Y,DTM.Z,RGB);
                axis equal;
                camlight
            end
            RM = makerefmat(min(X),min(Y),1, 1);
            W = refmatToWorldFileMatrix(RM);
            R = maprasterref(W,size(squeeze(RGB(:,:,1))),'postings');
            disp(filez(i).name);
            geotiffwrite(['/Volumes/FastDisk/RGB_AG/RGB/',filez(i).name(1:14),'RGB.tif'],RGB,R,'GeoKeyDirectoryTag', ...
                         info.GeoTIFFTags.GeoKeyDirectoryTag);
            %end
    end
end
return
% make subset for Fricktal
    cd /Volumes/lidarlab/data/KantonAargau/Mat_Files_Statistics/geotiffRGB
    filez = dir('*.tif');
    for i = 1:length(filez)
        nm = filez(i).name;
        east = str2num(nm(1:6));
        nort = str2num(nm(8:13));
        if east <= 654000 & nort >= 253000
            disp(['mv ',filez(i).name,' freakvalley']);
            unix(['mv ',filez(i).name,' freakvalley']);
        end
    end
    