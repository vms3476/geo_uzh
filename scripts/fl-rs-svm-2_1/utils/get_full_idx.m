function idx2=get_full_idx(indice)


temp=0;

app=[];
test=[];
for j=1:length(indice{1})
    app=[];
    test=[];
    temp=0;
    for i=1:length(indice)

        tapp=indice{i}(j).app ;
        ttest=indice{i}(j).test;
        app=[app;tapp+temp];
        test=[test;ttest+temp];
        temp=max(temp,max(app));
        temp=max(temp,max(test));


    end

    idx2(j).app=app;
    idx2(j).test=test;
end



