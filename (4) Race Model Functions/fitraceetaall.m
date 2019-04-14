function stable = fitraceetaall(data)

subjects = data(:,1);

s = max(subjects);

y = 1;

for c = 1:s
    if c~=3
    
    [fvalue,fit] = raceqeta(data,c)
    
    fvalue = fvalue;
    fit = fit;
    
    answer(y,1:7) = [c fvalue fit];
    
    colnames = {'Subject', 'BIC', 'mu1', 'mu2','sigma1', 'sigma2', 'eta'};
    
    stable = array2table(answer,'VariableNames',colnames);
    y= y+1;
end
end
    