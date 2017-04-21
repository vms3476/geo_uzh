function [xmin,vars,costmin,LOG]=nesterov_solver(x0,vars,data,opts)
% [x,cost,vars]=iterative_solver(x0,vars0,data,opts)
% function doing function executing iteratively an update function until
% convergence conditions
% input:
%   x0 : value of x before the descent (x is a column vector)
%   vars0: stucture needed by the cost function and the gradient function
%          (this structure can be updated for warm start or not )
%   data : constant data used in the cost and gradient function
%   opts : options fior the gradient descent
%   opts.q_func: function for the cost that is used like that :
%           [costnew,varsnew]=q_func(W,Wt,vars,data);
%   opts.minq_func: function for the cost that is used like that :
%           [x,costnew,varsnew]=min_func(x_1,vars,data);
%   opts.costfunc: function for the cost that is used like that :
%           [costnew,varsnew]=costfunc(x_1,vars,data);
%   opts.verbose: 0 if not, 1 if print + nbspaces
%   opts.log : 1 if saving LOG
%   opts.stopvarx : x variation threshold for stoping
%   opts.stopvarj : J cost variation threshold for stoping
%   opts.nbitermax: max number of iterations
%   opts.numericalprecision : numerical precision for the minimization
%
%  Author : Remi Flamary remi (dot) flamary (at) gmail (dot) com

if ~isfield(opts,'verbose'); opts.verbose=0; end
if ~isfield(opts,'log'); opts.log=1; end
if ~isfield(opts,'stopvarx'); opts.stopvarx=0; end
if ~isfield(opts,'stopvarj'); opts.stopvarj=1e-3; end
if ~isfield(opts,'nbitermax'); opts.nbitermax=100; end
if ~isfield(opts,'numericalprecision'); opts.numericalprecision=1e-5; end
if ~isfield(opts,'printfunc'); opts.printfunc=@print_0; end
if ~isfield(opts,'nu'); opts.nu=2; end
if ~isfield(opts,'L0'); opts.L0=1; end
if ~isfield(opts,'func_stop'); opts.func_stop=@func_0; end
if ~isfield(opts,'backtrack'); opts.backtrack=0; end


LOG=[];
J=[];
time=[];



W=x0;

% initialization
[J,vars]=opts.costfunc(W,vars,data);

if opts.verbose==1
%    opts.printfunc(0,opts.verbose,J,0,vars); 
%QUESTA FUNZIONE PRINTA TUTTO. DISATTIVATA X SCAZZO
end

if opts.log
    LOG.vars=vars;
end

nbverbose=0;
loop=1;
nbloop=1;
if ~isfield(opts,'L');
vars.L=norm(data.X'*data.X)*opts.C; % dont forget the regularization term
%Xt=data.X.*repmat(data.Y,1,size(data.X,2));
%vars.L=norm(Xt'*Xt)*opts.C; % dont forget the regularization term
else
vars.L=opts.L;    
end
vars.V=W;
vars.a_1=1;
vars.a=1;
vars.t=0;
vars.dt=zeros(size(W));

% MAIN LOOP

%tic;
W_1=W+1;

while loop


    % computing proximal operator on V
    if opts.log || opts.verbose|| opts.stopvarj~=0
        [W,ctemp,vars]=opts.minq_func(vars.V,vars,data);
        F=opts.costfunc(W,vars,data);
        J(nbloop+1)=F;
    end

    % if backtrack
    if opts.backtrack || loop==1
        while (F>opts.q_func(W,vars.V,vars,data))
            vars.L=vars.L*opts.nu;
            if opts.log || opts.verbose 
                [W,ctemp,vars]=opts.minq_func(vars.V,vars,data);
                F=opts.costfunc(W,vars,data);
            end
            J(nbloop+1)=F;
        end
    end

    %-------------------------------------------------
    % update
    %-----------------------------------------------
    %       Update ?? la Tseng
    %        vars.a=2/(vars.t-1+3);
    %
    %        vars.dt=W-W_1;
    %
    %        vars.V=W+(1-vars.a_1)/vars.a_1*vars.a*vars.dt;
    %
    %
    %        vars.a_1=vars.a;
    %        vars.t=vars.t+1;

    %-------------------------------------------------
    %      Real Fista
    %-------------------------------------------------
    vars.dt=W-W_1;

    vars.a_1=(1+sqrt(1+4*vars.a^2))/2;
    vars.V=W + (vars.a-1)/vars.a_1*vars.dt;
    vars.a=vars.a_1;



    if opts.verbose ~= 0
        % verbose print function
        %nbverbose=opts.printfunc(nbverbose,opts.verbose,J(nbloop+1),(J(nbloop+1)-J(nbloop))/J(nbloop),vars);
        %QUESTA FUNZIONE PRINTA TUTTO. DISATTIVATA X SCAZZO
    end


    if nbloop>opts.nbitermax
        loop=0;
        if opts.verbose ~= 0   disp('Max number of iteration reached'); end
    end


    if opts.func_stop(W,vars,data)
        loop=0;
        if opts.verbose ~= 0   disp('Stop condition reached'); end
    end

    if abs(J(nbloop+1)-J(nbloop)/J(nbloop))<=opts.stopvarj
        loop=0;
        if opts.verbose ~= 0     disp('DeltaJ convergence'); end
    end

    % x variation stopping criterion
    if max(max(max(abs(W-W_1))))/sum(sum(sum(abs(W))))<opts.stopvarx
        loop=0;
        if opts.verbose ~= 0     disp('Delta x convergence'); end
    end

    W_1=W;

    nbloop=nbloop+1;


end

xmin=W;
costmin=J(end);

if opts.log
    LOG.J=J;
    %        LOG.xtot=xtot;
    LOG.time=time;

end


end

%--------------------------------------------------------------------------
%           Sous Fonction
%--------------------------------------------------------------------------

function [nbverbose]=print_0(nbverbose,verbose,J,deltaJ,vars)

nbverbose = nbverbose+1;
if nbverbose == 20 || (nbverbose==1 && deltaJ==0);
    fprintf(1,[repmat('    ',1,verbose-1) '      Cost     Delta Cost \n']);
    nbverbose = 0;
end


fprintf(1,[repmat('    ',1,verbose-1) '| %11.4e | %8.4f |\n'],[J deltaJ ]);

end


function stop=func_0(W,vars,data)
    stop=0;
end