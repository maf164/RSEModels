function stable = racemodelmean(params,model)

y = 1;
for c = 1:14
    if(model=="race")
    mu1 = params(c, 2);
    mu1 = table2array(mu1);
    
    mu2 = params(c, 3);
    mu2 = table2array(mu2);
    
    sigma1 = params(c, 4);
    sigma1 = table2array(sigma1);
    
    sigma2 = params(c, 5);
    sigma2 = table2array(sigma2);
    
     eta = params(c, 6);
     eta = table2array(eta);
    
     rho = params(c, 7);
     rho = table2array(rho);
     
    elseif(model == "eta")
       mu1 = params(c, 2);
    mu1 = table2array(mu1);
    
    mu2 = params(c, 3);
    mu2 = table2array(mu2);
    
    sigma1 = params(c, 4);
    sigma1 = table2array(sigma1);
    
    sigma2 = params(c, 5);
    sigma2 = table2array(sigma2);
    
     eta = params(c, 6);
     eta = table2array(eta);
    
    elseif(model == "rho")
         mu1 = params(c, 2);
    mu1 = table2array(mu1);
    
    mu2 = params(c, 3);
    mu2 = table2array(mu2);
    
    sigma1 = params(c, 4);
    sigma1 = table2array(sigma1);
    
    sigma2 = params(c, 5);
    sigma2 = table2array(sigma2);
    
     rho = params(c, 6);
     rho = table2array(rho);
    elseif(model == "raab")
         mu1 = params(c, 2);
    mu1 = table2array(mu1);
    
    mu2 = params(c, 3);
    mu2 = table2array(mu2);
    
    sigma1 = params(c, 4);
    sigma1 = table2array(sigma1);
    
    sigma2 = params(c, 5);
    sigma2 = table2array(sigma2);
    
    end
    
  
x = 0:.001:1;


cdfa = laterCDF(x,mu1,sigma1);
idx = ~isnan(cdfa);
cdfa = cdfa(idx);
randval = rand(1,1000000);
[cdfa, mask] = unique(cdfa);
x = x(mask);
pa = interp1(cdfa,x,randval);
idx = ~isnan(pa);
pa = pa(idx);
meana = mean(pa);
sda = std(pa);

cdfv = laterCDF(x,mu2,sigma2);
idx = ~isnan(cdfv);
cdfv = cdfv(idx);
randval = rand(1,1000000);
[cdfv, mask] = unique(cdfv);
x = x(mask);
pv = interp1(cdfv,x,randval);
idx = ~isnan(pv);
pv = pv(idx);
meanv = mean(pv);
sdv = std(pv);

if model == "race"
cdfav = raceCDF(x,[mu1 mu2], [sigma1 sigma2],rho ,eta);
idx = ~isnan(cdfav);
cdfav = cdfav(idx);
randval = rand(1,1000000);
[cdfav, mask] = unique(cdfav);
x = x(mask);
pav = interp1(cdfav,x,randval);
idx = ~isnan(pav);
pav = pav(idx);
meanav = mean(pav);
sdav = std(pav);

elseif model == "eta"
    cdfav = raceCDF(x,[mu1 mu2], [sigma1 sigma2],0 ,eta);
idx = ~isnan(cdfav);
cdfav = cdfav(idx);
randval = rand(1,1000000);
[cdfav, mask] = unique(cdfav);
x = x(mask);
pav = interp1(cdfav,x,randval);
idx = ~isnan(pav);
pav = pav(idx);
meanav = mean(pav);
sdav = std(pav);

elseif model == "rho"
        cdfav = raceCDF(x,[mu1 mu2], [sigma1 sigma2],rho ,0);
idx = ~isnan(cdfav);
cdfav = cdfav(idx);
randval = rand(1,1000000);
[cdfav, mask] = unique(cdfav);
x = x(mask);
pav = interp1(cdfav,x,randval);
idx = ~isnan(pav);
pav = pav(idx);
meanav = mean(pav);
sdav = std(pav);

elseif model == "raab"
    cdfav = raceCDF(x,[mu1 mu2], [sigma1 sigma2],0 ,0);
idx = ~isnan(cdfav);
cdfav = cdfav(idx);
randval = rand(1,1000000);
[cdfav, mask] = unique(cdfav);
x = x(mask);
pav = interp1(cdfav,x,randval);
idx = ~isnan(pav);
pav = pav(idx);
meanav = mean(pav);
sdav = std(pav);
end

gain = getGain([sampleDown(transpose(pa),900000) sampleDown(transpose(pv),900000) sampleDown(transpose(pav),900000)]);

answer(y,1:7) = [meana sda meanv sdv meanav sdav gain];

stable = array2table(answer);
y = y+1;
end
end



