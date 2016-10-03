file=dir('*.mat');
% all.x=[];
% all.y=[];
% all.z=[];
% all.h=[];
% all.cbh=[];
% all.vol=[];
% all.dia=[];
% all.dia_EW=[];
% all.dia_NS=[];
% all.spec=[];
% 
% 
% for i=1:length(file)
%    load(file(i).name)
%    all.x=[all.x trees.x];
%    all.y=[all.y trees.y];
%    all.z=[all.z trees.z];
%    all.h=[all.h trees.h];
%    all.cbh=[all.cbh trees.cbh];
%    all.vol=[all.vol trees.vol];
%    all.dia=[all.dia trees.dia];
%    all.dia_EW=[all.dia_EW trees.dia_EW];
%    all.dia_NS=[all.dia_NS trees.dia_NS];
%    all.spec=[all.spec trees.spec];
%     
% end
for k=1:length(file)
    file_name=file(k).name;
    load(file_name)
    savename=file_name(1:13);
    
    [Einzelbaeume(1:length(trees.x)).Geometry]=deal('Point');
    for i=1:length(trees.x)
        
        Einzelbaeume(i).X=trees.x(i);
        Einzelbaeume(i).Y=trees.y(i);
        Einzelbaeume(i).ID=i-1;
        Einzelbaeume(i).Tree_z=trees.z(i);
        Einzelbaeume(i).Tree_h=trees.h(i);
        Einzelbaeume(i).cbh=trees.cbh(i);
        Einzelbaeume(i).Volume=trees.vol(i);
        Einzelbaeume(i).mean_dia=trees.dia(i);
        Einzelbaeume(i).dia_EW=trees.dia_EW(i);
        Einzelbaeume(i).dia_NS=trees.dia_NS(i);
        if trees.spec(i)==1;
            Einzelbaeume(i).Spec='Nadelbaum';
        else
            Einzelbaeume(i).Spec='Laubbaum';
        end
        
    end
    
    shapewrite(Einzelbaeume,strcat(savename,'_point'))
    
    [Baumumrisse(1:length(trees.x)).Geometry]=deal('Polygon');
    for j=1:length(trees.x)
        [Baumumrisse(j).X, Baumumrisse(j).Y]=outline(trees.olx{j}', trees.oly{j}');
        Baumumrisse(j).ID=j-1;
        if trees.spec(j)==1;
            Baumumrisse(j).Spec='Nadelbaum';
        else
            Baumumrisse(j).Spec='Laubbaum';
        end
        
    end
    shapewrite(Baumumrisse,strcat(savename,'_poly'))
end