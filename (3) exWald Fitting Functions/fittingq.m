function [fval,bestfit] = fittingq(data,model,s)

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


%Assign Quantiles
quantileA = quantile(A/1000,[0 .1,.3,.5,.7,.9]);
quantileV = quantile(V/1000,[0 .1,.3,.5,.7,.9]);
quantileAV = quantile(AV/1000,[0 .1,.3,.5,.7,.9]);

quants = [quantileA quantileV quantileAV];

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

%Fitting Procedure
fit = [];
options = optimset('MaxFunEvals',10000,'TolX',.00001,'PlotFcns',@optimplotfval);


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
function cost = minim(param,model)
    
emppropA = [ .1 .2 .2 .2 .2 .1];
emppropV = [ .1 .2 .2 .2 .2 .1];
emppropAV = [ .1 .2 .2 .2 .2 .1];

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
    
    modelpropA = modelquants(quants(1:6),ma,aa,ta);
    modelpropV = modelquants(quants(7:12),mv,av,tv);
    modelpropAV = modelquants(quants(13:18),mav,aav,tav);
    
   
logA = LogL(NA,emppropA,modelpropA,fpa);
logV = LogL(NV,emppropV,modelpropV,fpv);
logAV = LogL(NAV,emppropAV,modelpropAV,fpav);
totallog = (logA + logV + logAV);

[AIC, BIC] = aicbic(totallog,fpa+fpv+fpav,NA+NV+NAV);

cost = BIC;
end


%Calculate QMLE log-liklihood
    function logl = LogL(N,empiricalprops,modelprops,df)
        
        logl = sum(N.*empiricalprops.*log(modelprops));
    end

%For a paramater vector, evaulate the model quantiles at some rts -
%Zeheitleitner et al (2015)
    function res = modelquants(rts,m,a,t)
ex1 = exwaldcdf(rts,m,1,a,t);
ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
res = ex3;
    end
end



