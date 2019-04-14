function subject(data,c)

%Extracts the Visual, Auditory, and Audio-Visual for a single participant 
%as three seperate vectors from the full data matric. First input is data and
%should be  arranged in 3 columns with  subject number first, 
%RT second and condition third. Second input is participant number. 

%Written by Matthew Lansdell, University of St Andrews
%November 28 2018

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

assignin('base','A',A);
assignin('base','V',V);
assignin('base','AV',AV);