function stable = fitraceraaball(data)

subjects = data(:,1);

s = max(subjects);

y = 1;

for c = 1:s
    if c~=3
    
    [fvalue,fit] = raceqraab(data,c)
    
    fvalue = fvalue;
    fit = fit;
    
    answer(y,1:6) = [c fvalue fit];
    
    colnames = {'Subject', 'BIC', 'mu1', 'mu2','sigma1', 'sigma2'};
    
    stable = array2table(answer,'VariableNames',colnames);
    y= y+1;
end
end
    