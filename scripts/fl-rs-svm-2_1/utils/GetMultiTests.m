function res=GetMultiTests(opt)
% function returning a struct with the results from serveral file

res=struct( opt.varname , opt.var); 

for i=1:length(opt.var)

    optt=opt.opt;
    optt=setfield(optt,opt.varname,opt.var(i));
    optt=setfield(optt,'var',opt.var(i));
    filename=[opt.filename(optt) '.mat'];
    
    for j=1:length(opt.getvars)
        for k=1:j-1
        eval(['temp=load(''' filename ''',''' opt.getvars{j} ''');']);
        eval(['temp2=load(''' filename ''',''' opt.getvars{k} ''');']);
        
        for l=1:length(opt.tests)
        eval(['temp3=' opt.tests{l} '(temp.' opt.getvars{j} ',temp2.' opt.getvars{k} ');']);
        eval([' res.' opt.tests{l} '_' opt.getvars{j} '_' opt.getvars{k} '(' num2str(i) ')=temp3;']); 
        end
        end
    end
    
    

    
end