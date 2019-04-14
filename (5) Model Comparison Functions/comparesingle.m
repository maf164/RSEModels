function stable = comparesingle(data,condition)

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


A = A/1000;

rA = 1./A;

V = sort(V);


V = V/1000;

rV = 1./V;

if condition == 'A'
    
quantileA = quantile(A,[0 .1,.3,.5,.7,.9]);

quantilerA = quantile(rA,[0 .1,.3,.5,.7,.9]);

NA = size(A,1);
elseif condition == 'V'
    quantileA = quantile(V,[0 .1,.3,.5,.7,.9]);

quantilerA = quantile(rV,[0 .1,.3,.5,.7,.9]);

NA = size(V,1);
end

paramsr = [mean(rA) std(rA)];
options = optimset;
bestfit = inf;
for b = 1:10
try
[x,fval] = fminsearch(@(x) minimr(x),paramsr,options);

race=fval;
end
end
 for b=1:10
     try
ma = .01+ (30-.01).*rand(1,1);
aa = .01+ (20-.01).*rand(1,1);
ta = .01+ (5-.01).*rand(1,1);

paramsp = [ma aa ta];
    LB     = [0+eps 0+eps  0+eps ];  % These parameter values can't go below 0
UB     = [+Inf   +Inf   45  ];
[x,fval] = fminsearchbnd(@(x) minimp(x),paramsp,LB,UB,options);

if fval<bestfit;
    bestfit = fval;
pool=bestfit;
end
     end
 end
answer(y,1:3) = [c race pool];
    
colNames = {'Subject','Race','Pool'};
        stable = array2table(answer,'VariableNames',colNames);
        y = y+1;
        
 end
end


function costr = minimr(paramsr)

mu = paramsr(1);

sigma = paramsr(2);


empiricalprops = [ .1 .2 .2 .2 .2 .1];

ex1  = normcdf(quantilerA,mu,sigma);
ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
modelprops = ex3;

logr = sum(NA.*empiricalprops.*log(modelprops))

[AIC costr] = aicbic(logr,2,NA);
costr = costr;
end

function costp = minimp(paramsp)

ma = paramsp(1);

aa = paramsp(2);

ta = paramsp(3);

empiricalprops = [ .1 .2 .2 .2 .2 .1];

ex1 = exwaldcdf(quantileA,ma,1,aa,ta);
ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
res = ex3;
modelprops = res;

logp = sum(NA.*empiricalprops.*log(modelprops));


[AIC costp] = aicbic(logp,3,NA);
costp=costp;
end


end