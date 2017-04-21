function [acc,auc,svm,obj,LOG]=inflinearsvm_cv(xapp,yapp,xtest,ytest,indiceval,Cvec,nbCV,ratiocv,options,optionsFeatSel,InfoFeatures);


disp('CV phase')
    for itercv=1:nbCV
        %[xa,ya,xv,yv]=CreateDataAppVal(xapp,yapp,ratiocv);
        xa=xapp(indiceval(itercv).app,:);
        ya=yapp(indiceval(itercv).app);
        xv=xapp(indiceval(itercv).test,:);
        yv=yapp(indiceval(itercv).test); 
        for iC=1:length(Cvec)
            
            options.C=Cvec(iC);
            [svm,obj,LOG]=inflinearsvm(xa,ya,options,optionsFeatSel,InfoFeatures);
            
            ypred=inflinearsvmval(xv,svm);
            
            moyClassif=mean(sign(ypred)==yv);
            bccv(itercv,iC)=moyClassif;
            auccv(itercv,iC)=svmroccurve(ypred,yv);
        end;
    end;
    
    LOG.bccv=bccv;


    Mbccv=mean(bccv,1);
    [maxi,indC]=max(Mbccv);
    Copt=Cvec(indC);

    % learning phase
    disp('learning phase')
    options.C=Copt;
    
    tic
    [svm,obj,LOG2]=inflinearsvm(xapp,yapp,options,optionsFeatSel,InfoFeatures);
    LOG.timetot=toc;
    LOG.LOG2=LOG2;
    LOG.Copt=Copt;
    LOG.options=options;
    
    ypred=inflinearsvmval(xtest,svm);
    
    moyClassif=mean(sign(ypred)==ytest);
    acc=moyClassif;
    auc=svmroccurve(ypred,ytest);
    