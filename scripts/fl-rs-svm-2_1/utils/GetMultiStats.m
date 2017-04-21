function res=GetMultiStats(opt)
% function returning a struct with the stat results

res=struct( opt.varname , opt.var); 

if ~isfield(opt,'varsout')
    opt.varsout=opt.getvars;
end
if ~isfield(opt,'convert')
    opt.convert=0;
end
for i=1:length(opt.var)
    
    optt=opt.opt;
    optt=setfield(optt,opt.varname,opt.var(i));
    optt=setfield(optt,'var',opt.var(i));
    
    filename=[opt.filename(optt) '.mat'];
    
    for j=1:length(opt.getvars)
        eval(['temp=load(''' filename ''',''' opt.getvars{j} ''');']);
        if opt.convert
        for k=1:opt.nbpts
            eval(['temp2(k)=temp.' sprintf(opt.adress{j},k)  ';']);            
        end
        else
          eval(['temp2=temp.'  opt.getvars{j} ';']);     
        end
        
        for k=1:length(opt.funcs)
        eval([' res.' opt.funcs{k} '_' opt.varsout{j} '(' num2str(i) ')=' opt.funcs{k} '(temp2);']); 
        end
        
    end 
    
end