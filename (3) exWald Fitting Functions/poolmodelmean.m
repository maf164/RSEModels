function stable = poolmodelmean(params)

y = 1;
for c = 1:14
    
    ma = params(c, 2);
    ma = table2array(ma);
    
    aa = params(c, 3);
    aa = table2array(aa);
    
    ta = params(c, 4);
    ta = table2array(ta);
    
    mv = params(c, 5);
    mv = table2array(mv);
    
    av = params(c, 6);
    av = table2array(av);
    
    tv = params(c, 7);
    tv = table2array(tv);
    
    mav = params(c, 8);
    mav = table2array(mav);
    
    aav = params(c, 9);
    aav = table2array(aav);
    
    tav = params(c, 10);
    tav = table2array(tav);
   
x = 0:.001:1;

cdfa = exwaldcdf(x,ma,1,aa,ta);
randval = rand(1,1000000);
[cdfa, mask] = unique(cdfa);
x = x(mask);
pa = interp1(cdfa,x,randval);
idx = ~isnan(pa);
pa = pa(idx);
meana = mean(pa);
sda = std(pa);

cdfv = exwaldcdf(x,mv,1,av,tv);
randval = rand(1,1000000);
[cdfv, mask] = unique(cdfv);
x = x(mask);
pv = interp1(cdfv,x,randval);
idx = ~isnan(pv);
pv = pv(idx);
meanv = mean(pv);
sdv = std(pv);

cdfav = exwaldcdf(x,mav,1,aav,tav);
randval = rand(1,1000000);
[cdfav, mask] = unique(cdfav);
x = x(mask);
pav = interp1(cdfav,x,randval);
idx = ~isnan(pav);
pav = pav(idx);
meanav = mean(pav);
sdav = std(pav);

gain = getGain([sampleDown(transpose(pa),900000) sampleDown(transpose(pv),900000) sampleDown(transpose(pav),900000)]);

answer(y,1:7) = [meana sda meanv sdv meanav sdav gain];

stable = array2table(answer);
y = y+1;
end
end


