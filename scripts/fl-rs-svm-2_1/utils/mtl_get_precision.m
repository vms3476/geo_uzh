function [prec,log]= mtl_get_precision(Yte , Yval)
% USAGE function [prec_disc, prec_rej]= get_precision(vect_ytest,mat_ypred_val)
%
% vect_ytest est un vecteur donnant le label de chaque individu
% mat_ypred_val est la matrice donnant la confiance de l'appartenantce de chaque individu dans chaque classe



nbtasks=length(Yte);

for i=1:nbtasks
    
   log{i}= get_precision(Yte{i} , Yval{i});
   
   precs(i)=log{i}.prec;
    aucs(i)=log{i}.auc;
end
prec.log=log;
prec.precs=precs;
prec.aucs=aucs;
prec.prec=mean(precs);
prec.auc=mean(aucs);
prec.err=1-mean(precs);
%mat=matrice/sum(sum(matrice));
%precision=
