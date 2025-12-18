function ffagr=aggfactors
% program that computes aggregate factor returns from high-frequency data

% f: factors (T x k)
% h2: number of periods comprising a lower frequency

load numdaysquarterly.txt
num=numdaysquarterly;
%nl=length(num);
nl=181; % I replace the sum of all quarterly observations with a predetermined  number of 181 quarters. This number corresponds to 11391 daily observations.
load datareturns.txt
n=sum(num(1:nl));
y=datareturns(1:n,:);
k=5;
ff=y(:,1:k);
h=0;
for t=1:nl
    h2=num(t); %number of daily observations used to construct the lower-frequency measure
    hend=h+h2;
    ffagr(t,:)=sum(ff(h+1:hend,:)); % aggregate factors
    h=hend;
end
%mean(ffagr)
mean(ff)
std(ff)
%save aggfactors.txt ffagr -ascii
