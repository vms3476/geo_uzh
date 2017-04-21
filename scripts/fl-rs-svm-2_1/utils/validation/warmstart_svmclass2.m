function options=warmstart_svmclass2(X,options,model)

if ~isfield(model,'multiclass'); model.multiclass=0; end


if ~strcmp(options.kernel,'linear')

if model.multiclass==0
    options.alphainit=zeros(size(X,1),1);
    options.alphainit(model.pos)=abs(model.w);

else
    
    for i=1:model.nbclass
           options.alphainit{i}=zeros(size(X,1),1);
           options.alphainit{i}(model.models{i}.pos)=abs(model.models{i}.w);
        
    end
    
    
end

end