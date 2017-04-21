function res=GetMultiresults(opt)
% function returning a struct with the results from serveral file

if ~isfield(opt,'samenum'); opt.samenum=0; end
if ~isfield(opt,'twovars'); opt.twovars=0; end

res=struct( opt.varname , opt.var); 

res.var=opt.var;

for i=1:length(opt.var)
    
    optt=opt.opt;
    optt=setfield(optt,opt.varname,opt.var(i));
    if opt.samenum
      optt=setfield(optt,opt.varname2,opt.var(i));  
    end
    if opt.twovars
      optt=setfield(optt,opt.varname2,opt.var2(i));  
    end
   % optt=setfield(optt,'var',opt.var(i));
    filename=[opt.filename(optt) '.mat'];
    
    for j=1:length(opt.getvars)
        eval(['temp=load(''' filename ''',''' filterstr(opt.getvars{j}) ''');']);
        eval([' res.' opt.getvars{j} '(' num2str(i) ')=temp.' opt.getvars{j} ';']); 
        
    end 
    
end

end

function res=filterstr(str )

    k=strfind(str, '.');
    
    if isempty(k)
        res=str;
    else
       res=str(1:k(1)-1);
    end

end