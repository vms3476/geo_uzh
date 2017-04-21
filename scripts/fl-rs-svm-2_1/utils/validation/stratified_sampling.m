function [xapp, yapp, xtest ,ytest, OIapp, OItest ]= stratified_sampling(X, Y, kfold, melange_on, object_index)
%
% APPEL : [xapp,yapp,xtest,ytest]= stratified_sampling(X, Y, kfold, melange_on)
%
%  Echantillonage de dataset en 2 sous enssemble en respectant les proportions de classes
%
%      kfold =            ratio voulu entre le nombre d'echantillon en test et en app
%                                                  si kfold >1 :
%                                                       kfold = 5   =>  4/5 App, 1/5 Test
%                                                       kfold = 3   =>  2/3 App, 1/3 Test
%                                                  si kfold <1 :
%                                                       kfold = 1/3 =>  1/3 App, 2/3 Test
%                                                       kfold = 1/5 =>  1/5 App, 4/5 Test
%      melange_on =  booleen activant un melange des echantillons
%   Facultatif :
%      object_index = vecteur d'indice specifiant l'appartenance
%                                   des echantillons a des meta-objets
%                              Le sampling est alors fait de telle sorte que
%                                   les meta-objets ne soient pas sectionnes
%                                 (evite d'avoir 1 partie des meta-objets dans chaque dataset)
%
% B.LABBE  02-2010
%

if nargin<5
    object_index=[];
else
    melange_on=0;
end

if nargin<4
    melange_on=0;
end

selection_app=[];
selection_test=[];
labels=unique(Y);
nb_class=length(labels);

if isempty(object_index)
%--------------------------------------------------------------------------
    for i = 1:nb_class
        ind_i_app = find(Y==labels(i));
        nb_sample_i = length(ind_i_app);
        %
        if melange_on
            reordonne_app = randperm(nb_sample_i);
        else
            reordonne_app = 1:nb_sample_i;
        end
        %
        if kfold>1
                nb_select_app = floor(max(kfold-1,1)/kfold*nb_sample_i);
        else
                nb_select_app = floor(kfold*nb_sample_i);
        end
        %
        selection_app = [ selection_app ; ind_i_app(reordonne_app(  1:nb_select_app  ))];
        selection_test  = [ selection_test ; ind_i_app(reordonne_app(   nb_select_app+1  :  end  ))];
    end
%--------------------------------------------------------------------------
else
%--------------------------------------------------------------------------
    for i = 1:nb_class
        ind_i_app = find(Y==labels(i));
        nb_sample_i =length(ind_i_app);
        list_oi = unique(object_index(ind_i_app));
        nb_object_i = length(list_oi);
        %
        if melange_on
            reordonne_app = randperm(nb_object_i);
        else
            reordonne_app = 1:nb_object_i;
        end
        %
        if kfold>1
                nb_select_app = floor(max(kfold-1,1)/kfold*nb_object_i);
        else
                nb_select_app = floor(kfold*nb_object_i);
        end
        %
        for j = 1:nb_select_app
        selection_app = [ selection_app ;  find(object_index==list_oi(reordonne_app(j)  ))];
        end
        for j = nb_select_app+1 : nb_sample_i
        selection_test  = [ selection_test ; find(object_index==list_oi(reordonne_app(j)  ))];
        end
    end
%--------------------------------------------------------------------------
end

xapp=X(selection_app,:);
yapp=Y(selection_app,:);
xtest=X(selection_test,:);
ytest=Y(selection_test,:);

if nargout==6 & nargin==5
    OIapp=object_index(selection_app);
    OItest=object_index(selection_test);
end

