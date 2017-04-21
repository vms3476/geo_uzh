function [X,Y,Xval,Yval,Xte,Yte,info]=loadEPFL(opts,i)


if nargin<2
    i=0;
end

if ~isfield(opts,'props'), opts.props=[0.66 0 0.33]; end


opts.props=opts.props/sum(opts.props);

if isfield(opts,'nbpts')
nbpts=opts.nbpts;
idx=floor(cumsum([1 nbpts]));
else
    props=opts.props;
    nbpts=[0 0 0 0];
    idx=[];
end

%idx(1)=1;
      


for i=1:length(opts.tasks)
    
    t=opts.tasks(i);
    
    load(['S' num2str(t) '.mat'])
    
    Xin=x;
    Yin=y;
    
    
    if nbpts(3)==-1
        idx(4)=length(Yin);
    end
    
    info.nbpts(i)=length(Yin);
    
    if ~isfield(opts,'nbpts')
        nbmax=length(Yin);
         idx=floor(nbmax*cumsum([0 props]));
         idx(1)=1;
    end
    
X{i}=Xin(idx(1):idx(2)-1,:);
Y{i}=Yin(idx(1):idx(2)-1,:);
Xval{i}=Xin(idx(2):idx(3)-1,:);
Yval{i}=Yin(idx(2):idx(3)-1,:);
Xte{i}=Xin(idx(3):idx(4)-1,:);
Yte{i}=Yin(idx(3):idx(4)-1,:);

end

