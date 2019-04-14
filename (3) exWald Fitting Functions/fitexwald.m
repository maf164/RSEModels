function [fvalue,fit] = fitexwald(data,model,s)

%FITS TO QUANTILES


%Uses "fitting" to fit the ex-Wald model with several different starting
%parameters. Input arguments are data arranged in two columns with RT in column 1 and condition (A=1, V=2,
%AV=3) in column 2 and model which can be "drift", "driftter","ter",or
%"full".

%Written by Matthew Lansdell, University of St Andrews
%November 16 2018

%Index Data 
idx=data(:,1)<20;
    
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


cpA = getCP(A);
cpV = getCP(V);
cpAV = getCP(AV);
hold on
plot(A,cpA,'.')
plot(V,cpV,'.')
plot(AV,cpAV,'.')

%Default number of starting values to be tested
nStart = 10;
fit = [];
fvalue = Inf;

%Iterates through fitting function with random starting values and updates
%the fit if a lower BIC is obtained
for ss=1:nStart
    
try
[fval,bestfit] = fittingq(data,model,s);

fval = fval;
bestfit = bestfit;

if fval<fvalue
    fvalue = fval;
    fit = bestfit;
end

end
disp(fvalue)

end

%Output CDF graphs to compare fit 

Aparams = [bestfit(1) bestfit(2) bestfit(3)];
Vparams = [bestfit(4) bestfit(5) bestfit(6)];
AVparams = [bestfit(7) bestfit(8) bestfit(9)];

x = 0:.001:1;
xx= x*1000;

Acdf = exwaldcdf(x,Aparams(1),1,Aparams(2),Aparams(3));
Vcdf = exwaldcdf(x,Vparams(1),1,Vparams(2),Vparams(3));
AVcdf = exwaldcdf(x,AVparams(1),1,AVparams(2),AVparams(3));
hold on 
plot(xx,Acdf)
plot(xx,Vcdf)
plot(xx,AVcdf)

 legend('EmpiricalA','EmpiricalV','EmpiricalAV','ExWaldModelA','ExWaldModelV','ExWaldModelAV','location','southeast')
end
