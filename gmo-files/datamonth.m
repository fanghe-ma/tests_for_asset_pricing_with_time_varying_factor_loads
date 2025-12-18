function [days,aux]=datamonth(k)

load auxi.txt

n=length(auxi); j=0; days=0; freq=0;
for i=1:n
    j=j+1;    
    if (auxi(i)==1000)
        freq=freq+1;
        if (freq==k)
            days=days+1;
            aux(days)=j-k;
            j=0;
            freq=0;
        end
    end
end
days
save numdaysquarterly.txt aux -ascii
        