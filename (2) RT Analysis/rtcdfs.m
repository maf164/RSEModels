function rtcdfs(data)

%Analyzes all subject data at once and creates CDFs based on 6 quantiles.
%Data should be arranged in 3 columns with subject number first, RT second
%and condition third 

subjects = data(:,1);

s=max(subjects);

quantP = [.1,.3,.5,.7,.9];

for c=1:s
    
    idx=data(:,1)==c;
    
    x = data(idx,2:3);
    
    idx=x(:,2)==1;

A= x(idx);

idx=x(:,2)==2;

V= x(idx);

idx=x(:,2)==3;

AV= x(idx);

if not(isempty(AV))
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

AV = quantile(AV,quantP);
V = quantile(V,quantP);
A = quantile(A, quantP );

plot(A,quantP,'o')
plot(V,quantP,'o')
plot(AV,quantP,'o')

set(gca, 'xlim', [100 1000])

legend('Auditory','Visual','Audio-Visual','Location','southeast')

fname = sprintf('RTCDFS/rtcdfs%d.png', c);
saveas(gcf,fname)

    cla 
end
end