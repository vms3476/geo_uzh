function [acc,auc,svm,obj,LOG]=inflinearsvm_cv2(xapp,yapp,xtest,ytest,indiceval,Cvec,nbCSPvec,nbCV,ratiocv,options,optionsFeatSel,InfoFeatures);


disp('CV phase')
    for itercv=1:nbCV
        InfoFeaturesLoc=InfoFeatures;
        %[xa,ya,xv,yv]=CreateDataAppVal(xapp,yapp,ratiocv);
        xa=xapp(indiceval(itercv).app,:);
        ya=yapp(indiceval(itercv).app);
        xv=xapp(indiceval(itercv).test,:);
        yv=yapp(indiceval(itercv).test); 
        for iC=1:length(Cvec)
            for iCSP=1:length(nbCSPvec)
            
            InfoFeaturesLoc.CSPkeep=nbCSPvec(iCSP);
            options.C=Cvec(iC);
            [svm,obj,LOG]=inflinearsvm(xa,ya,options,optionsFeatSel,InfoFeaturesLoc);
            
            ypred=inflinearsvmval(xv,svm);
            
            moyClassif=mean(sign(ypred)==yv);
            bccv(itercv,iC,iCSP)=moyClassif;
            auccv(itercv,iC,iCSP)=svmroccurve(ypred,yv);
            end
        end;
    end;


    Mbccv=reshape(mean(bccv,1),size(bccv,2),size(bccv,3));
    maxi=max(max(Mbccv));
    [indC,indCSP]=find(Mbccv==maxi);
    Copt=Cvec(indC(1));
    CSPopt=nbCSPvec(indCSP(1));

    LOG.Copt=Copt;
    LOG.CSPopt=CSPopt;
    
    % learning phase
    disp('learning phase')
    options.C=Copt;
    InfoFeatures.CSPkeep=CSPopt;
    tic
    [svm,obj,LOG]=inflinearsvm(xapp,yapp,options,optionsFeatSel,InfoFeatures);
    LOG.timetot=toc;
    
    ypred=inflinearsvmval(xtest,svm);
    
    moyClassif=mean(sign(ypred)==ytest);
    acc=moyClassif;
    auc=svmroccurve(ypred,ytest);
    