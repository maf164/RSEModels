function [bestfval, fit] = raceqeta(data,s)

%This function fits the Race Model to emprical data through the use of
%Quantile Maximum Liklihood Estimation. Output is the BIC value and the
%best fitting paramaters. 

%Code written by Matthew Lansdell, University of St Andrews 2019

%Based on a variety of sources including Otto and Zehetlinier 



nStart = 100;

nGrid = 100;


%Prepare emprical data
idx=data(:,1)==s;
    
x = data(idx,2:3);
    
idx=x(:,2)==1;

A= x(idx);

idx=x(:,2)==2;

V= x(idx);

idx=x(:,2)==3;

AV= x(idx);


A = A/1000;
A = sort(A);
hold on 
A = 1./(A);


V = V/1000;
V = sort(V);
V = 1./(V);


AV = AV/1000;
AV = sort(AV);
AV = 1./(AV);


data  = [A(:);V(:);AV(:)];

idx = [     ones( length(A), 1 ); ...
        2 * ones( length(V), 1 ); ...
        3 * ones( length(AV), 1 ) ];

NA = size(A,1);
NV = size(V,1);
NAV = size(AV,1);

%Assign Quantiles
quantileA = quantile(A,[0 .1,.3,.5,.7,.9]);
quantileV = quantile(V,[0 .1,.3,.5,.7,.9]);
quantileAV = quantile(AV,[0 .1,.3,.5,.7,.9]);

quants = [quantileA quantileV quantileAV];




%Prepare initial paramater values 

startL = zeros(1,4);
for cc=1:2
    startL(cc)   = mean( data(idx==cc) );
    startL(cc+2) =  std( data(idx==cc) );
end

startR = getStartVal( data, idx, startL, nStart, nGrid );

 bestfval = Inf;
fit = [];

%Fitting Procedure

for ss=1:nStart
    try
param = [ startL, startR(ss,:) ];
param = [param(1) param(2) param(3) param(4) param(5)];
mu1 = param(1);
mu2 = param(2);
sigma1 = param(3);
sigma2 = param(4);
eta = param(5);





 options = optimset; 
    
    [x,fval] = fminsearch(@(x) minim(x),param, options);
    
    if fval<bestfval
        bestfval = fval;
        
        fit = x;
    end

    
    end
end

%Plot the cdfs of the best fitting paramaters
% modela = laterCDF(0:.001:1,fit(1),fit(3));
% modelv = laterCDF(0:.001:1,fit(2),fit(4));
% 
% modelav = raceCDF(0:.001:1,[fit(1) fit(2)],[fit(3);fit(4)],0,fit(5));
% 
% plot(0:.001:1,modela);
% plot(0:.001:1,modelv);
% plot(0:.001:1,modelav);


%Function being minimized
function cost = minim(param)
 mu1 = param(1);
mu2 = param(2);
sigma1 = param(3);
sigma2 = param(4);
eta = param(5);

emppropA = [ .1 .2 .2 .2 .2 .1];
emppropV = [ .1 .2 .2 .2 .2 .1];
emppropAV = [ .1 .2 .2 .2 .2 .1];

ex1  = normcdf(quants(1:6),mu1,sigma1);
ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
modelpropA = ex3;


ex1 = normcdf(quants(7:12),mu2,sigma2);

ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
modelpropV = ex3;


ex1 = maxNormCDF(quants(13:18), [mu1 mu2], [sigma1 sigma2] + eta, 0 );

ex2 = ex1;
ex2(1) = [];
ex2(6) = 1;
ex3 = ex2-ex1;
modelpropAV = ex3;

%modelpropAV = modelquants(quants(13:18),mu1,mu2,sigma1,sigma2,eta,rho);
    
   
costA = BICa(NA,emppropA,modelpropA);
costV = BICa(NV,emppropV,modelpropV);
costAV = BICa(NAV,emppropAV,modelpropAV);
log = costA + costV + costAV;
[AIC BIC] = aicbic(log,5,NA+NV+NAV);
cost = BIC;
end



function bic = BICa(N,empiricalprops,modelprops)
        
        bic = sum(N.*empiricalprops.*log(modelprops));
    end


function startR = getStartVal( data, idx, startL, nStart, nGrid )
    
% Range of start values for RHO and ETA
rhoVal = linspace( -0.95, 0.95, nGrid );
etaVal = linspace( -0.2,  0.5,  nGrid );
[ rhoVal, etaVal ] = meshgrid( rhoVal, etaVal );

% Grid search
logL = zeros( nGrid );
for ii=1:nGrid^2
    logL(ii) = sum( log( modelpdf( data, startL(1), startL(2), startL(3), startL(4), rhoVal(ii), etaVal(ii), idx ) ) );
end
logLsorted = sort( logL(:), 'descend' );

% Select start values based on grid search
startR = zeros(nStart,2);
for ss=1:nStart
    iStart = logL==logLsorted(ss);
    startR(ss,:) = [ rhoVal(iStart) etaVal(iStart) ];
end
end
function y = modelpdf( data, mu1, mu2, sigma1, sigma2, rho, eta, idx )
        
% Single signal conditions
y1 = normpdf( data(idx==1), mu1, sigma1 );
y2 = normpdf( data(idx==2), mu2, sigma2 );

% Redundant signals condition
yR = maxNormPDF( data(idx==3), [mu1 mu2], [sigma1 sigma2] + eta, rho );

% Collate likelihoods
y = [ y1; y2; yR];

end
end