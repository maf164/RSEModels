function stable = compare(data)

subjects = data(:,1);

s=max(subjects);
y = 1;

for c=1:s

if c ~= 3
 fullrace  = raceq(data,c);
rho = raceqrho(data,c);
eta = raceqeta(data,c);
raab = raceqraab(data,c);
        
        [tfvalue] = fitexwald(data,'ter',c);
        [dfvalue] = fitexwald(data,'drift',c);
        [ffvalue] = fitexwald(data,'full',c);
        [dtfvalue] = fitexwald(data,'driftter',c);

answer(y,1:9) = [c fullrace rho eta raab tfvalue dfvalue ffvalue dtfvalue];


colNames = {'Subject','fullrace','rho','eta','raab','ter','drift','full','driftter'};
        stable = array2table(answer,'VariableNames',colNames);
        y = y+1;
else end
end
end