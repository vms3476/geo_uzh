function [svm,obj]=svmclass2(X,Y,options,alphainit)
% USAGE [svm]=svmclass2(X,Y,options)
%
% Support vector machine for CLASSIFICATION
% This routine classify the training set with a support vector machine
% using quadratic programming algorithm selected by options.solver
%
% INPUT
%
% Training set
%      X  		: input data
%      Y  		: output data
% parameters
%		options.C		: Bound on the lagrangian multipliers
%		options.lambda		: Conditioning parameter for QP method
%       options.solver : svm solver used (and svm tyype) (see svmclass_chooser)
%                           0 : SVM  monqp.m (Active Constraints, Canu)
%                           1 : L2 SVM : (gradient descent, Chapelle)
%                           2 : SVM  (SimpleSVM, Looslie)
%                           3 : SVM svqp2 (SMO, Leon Bottou)
%                           4 : SVM/L2 SVM libsvm (SMO Chih-Jen Lin)
%                           5 : L2 SVM : monqp.m (Active Constraints, Canu)
%                           6 : LP SVM : (Linear programing Mangasarian 1998)
%
%		options.kernel		: kernel  type. classical kernel are
%
%		Name			parameters
%		'poly'		polynomial degree
%		'gaussian'	gaussian standard deviation
%
%		for more details see svmkernel
%
%		options.kerneloption : parameters of kernel
%
%		for more details see svmkernel
%
% 		options.verbose : display outputs (default value is 0: no display)
%
%      options.Span    : span matrix for semiparametric learning
%               This vector is sized Nbapp*Nbbasisfunction where
%               phi(i,j)= f_j(x(i));
%
%
%
%
%
% OUTPUT
%
% svm.xsup	coordinates of the Support Vector
% svm.w      weight
% svm.w0		bias
% svm.pos    position of Support Vector
% obj    Value of Objective function
%
%
% see also svmreg, svmkernel, svmval, svmclass


%	01/12/2009 R. Flamary




if ~isfield(options,'lambda'); options.lambda=1e-10; end
if ~isfield(options,'kernel'); options.kernel='gaussian'; end
if ~isfield(options,'kerneloption'); options.kerneloption=1; end
if ~isfield(options,'C'); options.C=1; end
if ~isfield(options,'solver'); options.solver=0; end
if ~isfield(options,'span'); options.span=[]; end
if ~isfield(options,'verbose'); options.verbose=0; end
if ~isfield(options,'alphainit'); options.alphainit=[]; end
if ~isfield(options,'prob'); options.prob=0; end
if ~isfield(options,'sel'); options.sel=[]; end

if ~isfield(options,'multiclass'); options.multiclass=0; end

if length(unique(Y))>2
    options.multiclass=1;
end

if options.multiclass==0

    if nargin>3
        options.alphainit=alphainit;
    end

    options.C=options.C/length(Y);

    if strcmp(options.kernel,'linear')
        
        switch options.prob

            case 1 % lda

                [w,w0]=ldaclass(X,Y,1/options.C);
                obj=0;

            otherwise % svm
                [w,w0,obj]=lsvmclass(X,Y,options.C);

        end
              
        svm.xsup=[];
        svm.w=w;
        svm.w0=w0;
        svm.pos=[];
        svm.kerneloption=options.kerneloption;
        svm.kernel=options.kernel;
        svm.span=options.span;
        
  


    else
        [xsup,w,d,pos,timeps,alpha,obj]=svmclass_chooser(X,Y,options.C,options.lambda,options.kernel,options.kerneloption,options.verbose,options.span, options.alphainit,options.solver,options);


        svm.xsup=xsup;
        svm.w=w;
        svm.w0=d;
        svm.pos=pos;
        svm.kerneloption=options.kerneloption;
        svm.kernel=options.kernel;
        svm.span=options.span;
    end


else

    vals=sort(unique(Y),'ascend');

    if options.verbose~=0
        disp(['Multiclass SVM classification (' num2str(length(vals)) 'classes)' ]);
    end

    switch options.multiclass

        case 1 % one against all


            options2=options;
            options2.multiclass=0;

            for i=1:length(vals)

                if ~isempty(options.alphainit)
                    options2.alphainit=options.alphainit{i};
                end


                if options.verbose~=0
                    disp(['Classification :' num2str(i) '/' num2str(length(vals)) ]);
                end;
                Yt=(Y==vals(i))-(Y~=vals(i));
                [models{i},obj(i)]=svmclass2(X,Yt,options2);

            end

            svm.multiclass=options.multiclass;
            svm.nbclass=length(vals);
            svm.models=models;
            svm.kerneloption=options.kerneloption;
            svm.kernel=options.kernel;
            svm.vals=vals;




    end

end


function [xsup,w,d,pos,timeps,alpha,obj]=svmclass_chooser(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit,solver,options)
% USAGE [xsup,w,b,pos,timeps,alpha,obj]=svmclass(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit)
%
% Support vector machine for CLASSIFICATION
% This routine classify the training set with a support vector machine
% using quadratic programming algorithm (active constraints method)
%
% INPUT
%
% Training set
%      x  		: input data 
%      y  		: output data
% parameters
%		c		: Bound on the lagrangian multipliers     
%		lambda		: Conditioning parameter for QP method
%		kernel		: kernel  type. classical kernel are
%
%		Name			parameters
%		'poly'		polynomial degree
%		'gaussian'	gaussian standard deviation
%
%		for more details see svmkernel
% 
%		kerneloption : parameters of kernel
%
%		for more details see svmkernel
%
% 		verbose : display outputs (default value is 0: no display)
%
%     Span    : span matrix for semiparametric learning 
%               This vector is sized Nbapp*Nbbasisfunction where
%               phi(i,j)= f_j(x(i));
%    
%       solver: solver used : 0 : Active constraint (SVM-KM, Canu)
%                             1 : gradient descent (Chapelle)
%                             2 : SimpleSVM (Loosli)
%                             3 : svqp2 SMO (Bottou)
%
%
% OUTPUT
%
% xsup	coordinates of the Support Vector
% w      weight
% b		bias
% pos    position of Support Vector
% timeps time for processing the scalar product
% alpha  Lagragian multiplier
% obj    Value of Objective function
%
%
% see also svmreg, svmkernel, svmval

%   24/11/09 R. Flamary
%	21/09/97 S. Canu
%	04/06/00 A. Rakotomamonjy   -inclusion of other kernel functions
%	04/05/01 S. Canu            -inclusion of multi-constraint optimization for frame-SVM
%
%       scanu@insa-rouen.fr, alain.rakoto@insa-rouen.fr

if nargin <10
    solver=0;
end

if nargin <11
    options=0;
end

switch solver

    case 1
         [xsup,w,d,pos,timeps,alpha,obj]=svmclasschapelle(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);
         
     
    case 2
         [xsup,w,d,pos,timeps,alpha,obj]=svmclasssimple(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);
      
    case 3
         [xsup,w,d,pos,timeps,alpha,obj]=svmclass_bottou(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);
    
              
    case 4
         [xsup,w,d,pos,timeps,alpha,obj]=svmclasslib(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);

    case 5
        
         [xsup,w,d,pos,timeps,alpha,obj]=svmclassL2(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);

    case 6 
        
         [xsup,w,d,pos]=LPsvmclass(x,y,c,lambda,kernel,kerneloption,verbose);
        timeps=0;
        alpha=abs(w);
        obj=0;
        
    case 7
         [xsup,w,d,pos,timeps,alpha,obj]=svmclass_apg(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit,options);

        
    otherwise
         [xsup,w,d,pos,timeps,alpha,obj]=svmclass(x,y,c,lambda,kernel,kerneloption,verbose,span, alphainit);
     
         
end










