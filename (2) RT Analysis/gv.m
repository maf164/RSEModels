function stable = gv(data)

%Computes the RSE gain and RMI violations for all participants in a
%dataset. Also returns the means of each condition and a simpler method of
%calculating the RSE gain %Data should be arranged in 3 columns with 
%subject number first, RT second and condition third.


%Pulls out each subject data using logical indexing
subjects = data(:,1);

s=max(subjects);

answer = [];

y = 1;

for c=15:15
    
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
    
  %Ensures there are an equal number of observations for each conditions 
  %Downsamples if not
    if length(AV)~= length(V) |  length(AV)~= length(A) | length(V)~= length(A);
    
    l = [length(AV) length(V) length(A)];
   low =  min(l);
     
         
    A = sampleDown(A,low);
    V = sampleDown(V,low);
    AV = sampleDown(AV,low);
    end

         m = [A V AV];
%Uses RSEbox functions to calculate
         G = getGain(m);
         Vi = getViolation(m);
         ma = mean(A);
         mv = mean(V);
         mav = mean(AV);
         
         %Method used to calculate RSE as described by Zehetleinter et al
         %(2015)
         Z = min(ma,mv)-mav;

         answer(y,1:7) = [c G Vi ma mv mav Z];


         colNames = {'Subject','Gain','Violation','AMean','VMean','AVMean','Z'};
         stable = array2table(answer,'VariableNames',colNames);
         y = y+1;
    end
    
end

end