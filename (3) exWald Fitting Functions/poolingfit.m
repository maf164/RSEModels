function [fval,bestfit] = poolingfit(data,model,s)

%FITS TO QUANTILES

%Function that finds the best fitting paramaters for four different types
%of the ex-Wald model as specified by Zehetleitner et al(2015). Input
%arguments are data arranged in two columns with RT in column 1 and condition (A=1, V=2,
%AV=3) in column 2, and model which can be "drift", "driftter","ter",or
%"full". The different models change which paramaters are held constant
%across the three conditions. 

%Written by Matthew Lansdell, University of St Andrews
%November 16 2018

%Code is based on fitting procedure implemented by Zehetleitner et al
%(2015)


%Creates an index to pull out the 3 conditions from the data.
idx=data(:,1)==s;
    
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

A = A/1000;
V = V/1000;
AV = AV/1000;




NA = size(A,1);
NV = size(V,1);
NAV = size(AV,1);

%Functions for creating random initial starting paramaters
if(model=="drift")
    
    samplefn = samplepardrift;
    
elseif(model=="ter")
    
    samplefn = sampleparter;
    
elseif(model=="driftter")
    
    samplefn = samplepardriftter;
    
elseif(model=="full")
    
    samplefn = sampleparfull;
end

function y=sampleparfull

ma = .01+ (30-.01).*rand(1,1);
aa = .01+ (20-.01).*rand(1,1);
ta = .01+ (5-.01).*rand(1,1);
mv = .01+ (30-.01).*rand(1,1);
av = .01+ (20-.01).*rand(1,1);
tv = .01+ (5-.01).*rand(1,1);
mav = .01+ (30-.01).*rand(1,1);
aav = .01+ (20-.01).*rand(1,1);
tav = .01+ (5-.01).*rand(1,1);
y = [ma aa ta mv av tv mav aav tav];
end


function y = samplepardrift
ma = .01+ (30-.01).*rand(1,1);
mv = .01+ (30-.01).*rand(1,1);
mav = .01+ (30-.01).*rand(1,1);
a = .01+ (20-.01).*rand(1,1);
t = .01+ (5-.01).*rand(1,1);
y = [ma a t mv a t mav a t];
end

function y=sampleparter
m = .01+ (30-.01).*rand(1,1);
a = .01+ (20-.01).*rand(1,1);
ta = .01+ (5-.01).*rand(1,1);
tv = .01+ (5-.01).*rand(1,1);
tav = .01+ (5-.01).*rand(1,1);
y = [m a ta m a tv m a tav];

 
end

function y=samplepardriftter

ma = .01+ (30-.01).*rand(1,1);
a = .01+ (20-.01).*rand(1,1);
ta = .01+ (5-.01).*rand(1,1);
mv = .01+ (30-.01).*rand(1,1);
tv = .01+ (5-.01).*rand(1,1);
mav = .01+ (30-.01).*rand(1,1);
tav = .01+ (5-.01).*rand(1,1);
y = [ma a ta mv a tv mav a tav];
end

 %samplefn = [23.35	4.9	0.12	46.34	10.92	0.06	27.02	5.56	0.03
%];
%Fitting Procedure
fit = [];
options = optimset('MaxFunEvals',10000,'TolX',.00001,'PlotFcns',@optimplotfval);

%options = optimset('fminsearch','MaxFunEvals',100000,'MaxIter',100000,'TolFun',.00001,'TolX',.00001,'Display','iter')


if(model=="drift")
    
LB     = [0+eps 0+eps  0+eps 0+eps 0+eps ];  % These parameter values can't go below 0
UB     = [+Inf   +Inf   45  +Inf   +Inf ]; % These are heuristic starting values based on Palmer et al(2012)
    samplefn = [samplefn(1) samplefn(2) samplefn(3) samplefn(4) samplefn(7)] ;
    [x,fval,exitflag,output] = fminsearchbnd(@(x) minim(x,model),samplefn,LB,UB,options)
    
    fval = fval;
    bestfit = [x(1) x(2) x(3) x(4) x(2) x(3) x(5) x(2) x(3)];
    
elseif(model=="ter")
    LB     = [0+eps 0+eps  0+eps 0+eps 0+eps  ];  % These parameter values can't go below 0
UB     = [+Inf   +Inf   45   45     45 ];  % These are heuristic starting values based on Palmer et al(2012)
    samplefn = [samplefn(1) samplefn(2) samplefn(3) samplefn(6) samplefn(9)] ;
    [x,fval,exitflag,output] = fminsearchbnd(@(x) minim(x,model),samplefn,LB,UB,options)
    bestfit = [x(1) x(2) x(3) x(1) x(2) x(4) x(1) x(2) x(5)];
   
elseif(model=="driftter")
    LB     = [0+eps 0+eps  0+eps 0+eps 0+eps  0+eps 0+eps ];  % These parameter values can't go below 0
UB     = [+Inf   +Inf   45  +Inf     45  +Inf    45 ];  % These are heuristic starting values based on Palmer et al(2012)
    samplefn = [samplefn(1) samplefn(2) samplefn(3) samplefn(4) samplefn(6) samplefn(7) samplefn(9)] ;
    [x,fval,exitflag,output] = fminsearchbnd(@(x) minim(x,model),samplefn,LB,UB,options)
    bestfit = [x(1) x(2) x(3) x(4) x(2) x(5) x(6) x(2) x(7)];
    
    
elseif(model=="full")
    LB     = [0+eps 0+eps  0+eps 0+eps 0+eps  0+eps 0+eps 0+eps  0+eps];  % These parameter values can't go below 0
UB     = [+Inf   +Inf   45  +Inf   +Inf   45  +Inf   +Inf   45 ];  % These are heuristic starting values based on Palmer et al(2012)
    samplefn = samplefn;
    %samplefn = [23.35 4.9 .12 46.34 10.92 .06 27.02 5.56 .03]
    [x,fval,exitflag,output] = fminsearchbnd(@(x) minim(x,model),samplefn,LB,UB,options)
    bestfit = [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9)];
    
    
    
end    



%Function being minimized
function totallog = minim(param,model)
    

if(model=="drift")
    ma = param(1) ;
    aa = param(2);
    ta = param(3);
    mv = param(4);
    av = param(2);
    tv = param(3);
    mav = param(5);
    aav = param(2);
    tav = param(3);
    
    fpa = 3;
    fpv = 1;
    fpav =1;
elseif(model=="ter");
    ma = param(1) ;
    aa = param(2);
    ta = param(3);
    mv = param(1);
    av = param(2);
    tv = param(4);
    mav = param(1);
    aav = param(2);
    tav = param(5);
    
     fpa = 3;
    fpv = 1;
    fpav =1;
elseif(model=="driftter");
    ma = param(1) ;
    aa = param(2);
    ta = param(3);
    mv = param(4);
    av = param(2);
    tv = param(5);
    mav = param(6);
    aav = param(2);
    tav = param(7);
    
    
     fpa = 3;
    fpv = 2;
    fpav = 2;
elseif(model=="full")
    ma = param(1) ;
    aa = param(2);
    ta = param(3);
    mv = param(4);
    av = param(5);
    tv = param(6);
    mav = param(7);
    aav = param(8);
    tav = param(9);
     fpa = 3;
    fpv = 3;
    fpav = 3;
end
    % Number of interquantile intervals
m = 7; 

% Vector of cumulative probabilities 
p = getCP(m-1); 

% Vector of quantile estimates
qa = quantile(A, [0 .1,.3,.5,.7,.9]); 
qa1 = qa';
qv = quantile(V, [0 .1,.3,.5,.7,.9]); 
qv1 = qv';
qav = quantile(AV, [0 .1,.3,.5,.7,.9]); 
qav1 = qav';

% Vector of counts 
Na = histcounts(A, [-Inf; qa1;  Inf]); 

Nv = histcounts(V, [-Inf; qv1;  Inf]); 

Nav = histcounts(AV, [-Inf; qav1 ; Inf]); 


% LATER CDF
pa = [0 exwaldcdf(qa, ma,1,aa,ta) 1]; 
pv = [0 exwaldcdf(qv, mv,1,av,tv) 1]; 
pav = [0 exwaldcdf(qav, mav,1,aav,tav) 1]; 

% log likelihood
logLa = -sum(log((pa(2:m+1)-pa(1:m)) .^Na)); 
logLv = -sum(log((pv(2:m+1)-pv(1:m)) .^Nv)); 
logLav = -sum(log((pav(2:m+1)-pav(1:m)) .^Nav)); 

totallog = logLa +  logLv +  logLav;
end
end

  




