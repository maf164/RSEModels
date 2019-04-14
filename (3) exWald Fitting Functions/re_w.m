function z=re_w(x,y);
% Computes the real part of w(x+iy) for the ex-Wald.
%
% For use with EXWALDPDF
% 
 
% Adapted from the equations presented in:
%   Schwarz, W. (2001). The ex-Wald distribution as a descriptive model of
%   response time. Behavior Research Methods, Instruments, & Computers,
%   33(4), 457-469.
%
% Original Code from:
%   Wolf Schwarz
%   Dept. of Psychology
%   University of Potsdam
%   wschwarz@rz.uni-potsdam.de
%
% Adapted for MATLAB by: 
%   Evan McHughes Palmer
%   Visual Attention Laboratory
%   Harvard Medical School
%   palmer@search.bwh.harvard.edu
%   
%   September 21, 2005

% Constants for approximation to re[w] for x > 3.9 or y > 3.0
wa=[0.4613135 0.09999216 0.002883894];
wb=[0.1901635 1.7844927 5.5253437];

hi=find(x>3.9|y>3);
lo=find(x<=3.9&y<=3);

u=zeros(1,length(x));
v=zeros(1,length(x));
z=zeros(1,length(x));

x2=zeros(1,length(x));
y2=zeros(1,length(x));

x2=y;
y2=x;

fs=zeros(1,length(x));
fc=zeros(1,length(x));
f1=zeros(1,length(x));
f2=zeros(1,length(x));
zn=zeros(1,length(x));
zp=zeros(1,length(x));
ffn=zeros(1,length(x));
ggn=zeros(1,length(x));

%if (x>3.9|y>3)

    for i=1:3
        aa=wa(i);
        bb=wb(i);
        
        u(hi)=u(hi)+aa.*(x(hi).^2-y(hi).^2-bb)./((x(hi).^2-y(hi).^2-bb).^2+4.*x(hi).^2.*y(hi).^2);
        v(hi)=v(hi)-aa.*2.*x(hi).*y(hi)./((x(hi).^2-y(hi).^2-bb).^2+4.*x(hi).^2.*y(hi).^2);
    end;
    
    z(hi)=-x(hi).*v(hi)-y(hi).*u(hi);
    
%else
    
    og_n=20;    % Sets upper series summation limit
    
    % First, compute erf(y+ix) following 7.1.29 in Abramowitz & Stegun (1965)
    
    % Swap x and y
    % DID THIS ABOVE
    
    fs(lo)=sin(2.*x2(lo).*y2(lo));
    fc(lo)=cos(2.*x2(lo).*y2(lo));
    
    u(lo)=2.*normcdf(x2(lo).*sqrt(2),0,1)-1;
    u(lo)=u(lo)+(1-fc(lo)).*(exp(-x2(lo).*x2(lo)))./(2.*pi.*x2(lo));
    
    v(lo)=fs(lo).*(exp(-x2(lo).*x2(lo)))./(2.*pi.*x2(lo));
    
    f1(lo)=(2./pi).*exp(-x2(lo).*x2(lo));
    
    for n=1:og_n
        f2(lo)=1./(n.^2+4.*x2(lo).^2);
        
        zn(lo)=exp(-0.25.*n.*(n-4.*y2(lo)));
        zp(lo)=exp(-0.25.*n.*(n+4.*y2(lo)));
        
        ffn(lo)=2.*x2(lo).*exp(-0.25.*n.^2);
        ffn(lo)=ffn(lo)-x2(lo).*(zn(lo)+zp(lo)).*fc(lo);
        ffn(lo)=ffn(lo)+n.*0.5.*(zn(lo)-zp(lo)).*fs(lo);
    
        ggn(lo)=x2(lo).*(zn(lo)+zp(lo)).*fs(lo);
        ggn(lo)=ggn(lo)+n.*0.5.*(zn(lo)-zp(lo)).*fc(lo);
        
        u(lo)=u(lo)+f1(lo).*f2(lo).*ffn(lo);
        v(lo)=v(lo)+f1(lo).*f2(lo).*ggn(lo);
    end
    
    % Swap x and y
    x3=y2;
    y3=x2;
    
    % Next, compute w(x+iy) from erf(y+ix)
    
    fs(lo)=sin(2.*x3(lo).*y3(lo));
    fc(lo)=cos(2.*x3(lo).*y3(lo));
    
    z(lo)=(fc(lo).*(1-u(lo))+v(lo).*fs(lo)).*exp(y3(lo).^2-x3(lo).^2);
%end
 
 