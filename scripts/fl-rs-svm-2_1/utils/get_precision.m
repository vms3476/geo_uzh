function [prec]= get_precision(Yte , Yval)
% USAGE function [prec_disc, prec_rej]= get_precision(vect_ytest,mat_ypred_val)
%
% vect_ytest est un vecteur donnant le label de chaque individu
% mat_ypred_val est la matrice donnant la confiance de l'appartenantce de chaque individu dans chaque classe

if iscell(Yte)
    prec=mtl_get_precision(Yte , Yval);

else

    if size(Yval,2)>1

        %
        [temp,Ytemp]=max(Yval,[],2);
        prec=assessment(Yte,Ytemp,'class');
        prec2=get_precision_mc(Yte , Yval);
        prec.err=prec2.EC;
        prec.prec=prec2.CC;
        prec.pprec=mean(prec2.cc);
        prec.auc= prec.pprec;


        prec.mse=mean((Ytemp-Yte).^2);
        % prec.mse=

    else

        %labels=unique(vect_ytest);
        nb_indiv = length(Yte);
        %nb_class = length(labels);

        Ypred=(Yval>0)-(Yval<0);

        %---------------------------------------------
        %[temp, ycalc]=max(mat_ypred_val,[],2);
        % ycalc=labels(ycalc);
        %---------------------------------------------
        prec.prec=sum(Yte==Ypred)/nb_indiv;
        prec.err=1-prec.prec;
        %prec_disc.labels=labels;
        % for i=1:nb_class
        %   prec_disc.cc{i}=sum(ycalc==labels(i) & vect_ytest==labels(i))/ sum(vect_ytest==labels(i));
        % end
        labels=[-1 1];
        for i=1:2
            for j=1:2
                prec_disc.matrice(i,j)=sum(Yte==labels(i) & Ypred==labels(j));
            end
        end
        %---------------------------------------------
        [prec.auc,prec.tpr,prec.fpr]=svmroccurve(Yval,Yte);

    end
    %mat=matrice/sum(sum(matrice));
    %precision=

end
