function stable = rsegain(data)


subjects = data(:,1);

s = max(subjects);


y = 1;

for c=1:s
    
 if c ~=3


idx=data(:,1)==c;
    
x = data(idx,2:3);
    
idx=x(:,2)==1;

A= x(idx);

idx=x(:,2)==2;

V= x(idx);

idx=x(:,2)==3;

AV= x(idx);

A = sort(A);

V = sort(V);

AV = sort(AV);

A = sampleDown(A,200);
V = sampleDown(V,200);
AV = sampleDown(AV,200);

cpa = getCP(A);
cpv = getCP(V);
cpav = getCP(AV);

gricedata = [A V];
grice = getGrice(gricedata);
gricecp = getCP(grice);

gaindata = [A V AV];
gain = getGain(gaindata);
raab = getRaab([A V])
raab = getGain([A V raab])

answer(y,1:3) = [c gain raab];

colNames = {'Subject','Gain','Raab'};
        stable = array2table(answer,'VariableNames',colNames);
        y = y+1;
 else
 end
end
end





