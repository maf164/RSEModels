function rseviolationgraph(data)


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
hold on
plot(A,cpa,'.')
plot(V,cpv,'.')
plot(AV,cpav,'.')

gricedata = [A V];
grice = getGrice(gricedata);

miller = getMiller(gricedata);
millercp = getCP(miller);

plot(miller,millercp,'.')
 fill = [AV miller];
fillArea(fill,[ 0.8 0.3 0.1 ],1)


fname = sprintf('RSEviolation/rseviolation%d.png', c);
saveas(gcf,fname)

    cla
 else
 end
end
end







