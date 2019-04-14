function stable = fitexwaldall(data,model)

subjects = data(:,1);

s = max(subjects);

y = 1;

for c = 1:s
    if c~=3
    
    [fvalue,fit] = fitexwald(data,model,c)
    
    fvalue = fvalue;
    fit = fit;
    
    answer(y,1:11) = [c fvalue fit];
    
    colnames = {'Subject', 'BIC', 'ma', 'aa','ta', 'mv', 'av', 'tv', 'mav', 'aav', 'tav'};
    
    stable = array2table(answer,'VariableNames',colnames);
    y= y+1
end
end
    