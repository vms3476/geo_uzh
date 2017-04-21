function [Xapp, Yapp, Xvalid, Yvalid] = get_fold_k_crossVal(X, Y, propAppValid, offsetProp)
%
% APPEL : [Xapp, Yapp, Xvalid, Yvalid] = get_fold_k_crossVal(X, Y, propAppValid, offsetProp)
%
%
% B.LABBE 2 04 2010

% decalage des donnees sans recouvrement pour obtenir le K-ieme fold de la crossVal
melange_on=0;

[xtemp1, ytemp1, xtemp2, ytemp2]=stratified_sampling(X,Y,offsetProp, melange_on);

X=[xtemp2 ; xtemp1];
Y=[ytemp2 ; ytemp1];

% decoupage selon les proportions AppValid imposees
[Xapp, Yapp, Xvalid, Yvalid]=stratified_sampling(X,Y, propAppValid, melange_on);

