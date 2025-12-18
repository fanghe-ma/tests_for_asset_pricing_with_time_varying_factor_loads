function [rag,omega]=realizedmeasures
% program that computes realized covariance measures from high-frequency data
% r: cross-section of returns (T x 1)
% f: factors (T x k)
% h2: number of periods comprising a lower frequency

%load numdaysmonthly.txt
%num=numdaysmonthly;
load numdaysquarterly.txt
num=numdaysquarterly;
%sum(num(1:181))
nl=length(num);
%nl=181; % I replace the sum of all quarterly observations with a predetermined  number of 181 quarters. This number corresponds to 11391 daily observations.
load datareturns.txt
%n=length(datareturns);
n=sum(num(1:nl));
y=datareturns(1:n,:);
k=5;
m=47;
ff=y(:,1:k);

%plot(datareturns(1:11400,6))
%find(datareturns(:,6)<0.001)
%pause
% estimate of realized covariance
for i=1:m
    clear r; clear rag; 
    omega=zeros(nl,k); h=0;
    r=y(:,6+i)-y(:,6); % excess returns 
    rp=y(:,6+i)./y(:,6)-1; % percentage excess return
    for t=1:nl
        h2=num(t); %number of daily observations used to construct the lower-frequency measure
        hend=h+h2;
        omega(t,1:k)=r(h+1:hend)'*ff(h+1:hend,:); % (h x 1 )' x (h x k)    
        rag(t)=sum(r(h+1:hend)); % (1 x 1)
        ragp(t)=sum(rp(h+1:hend)); % (1 x 1)
        h=hend;
    end
    rag=rag'; ragp=ragp';
    ret(:,i)=rag;
    retp(:,i)=ragp;
    
    if (i==1)
        save omegareg1.txt omega -ascii
    elseif (i==2)
        save omegareg2.txt omega -ascii
    elseif (i==3)
        save omegareg3.txt omega -ascii
    elseif (i==4)
        save omegareg4.txt omega -ascii
    elseif (i==5)
        save omegareg5.txt omega -ascii
    elseif (i==6)
        save omegareg6.txt omega -ascii
    elseif (i==7)
        save omegareg7.txt omega -ascii
    elseif (i==8)
        save omegareg8.txt omega -ascii
    elseif (i==9)
        save omegareg9.txt omega -ascii
    elseif (i==10)
        save omegareg10.txt omega -ascii
    elseif (i==11)
        save omegareg11.txt omega -ascii
    elseif (i==12)
        save omegareg12.txt omega -ascii
    elseif (i==13)
        save omegareg13.txt omega -ascii
    elseif (i==14)
        save omegareg14.txt omega -ascii
    elseif (i==15)
        save omegareg15.txt omega -ascii
    elseif (i==16)
        save omegareg16.txt omega -ascii
    elseif (i==17)
        save omegareg17.txt omega -ascii
    elseif (i==18)
        save omegareg18.txt omega -ascii
    elseif (i==19)
        save omegareg19.txt omega -ascii
    elseif (i==20)
        save omegareg20.txt omega -ascii
    elseif (i==21)
        save omegareg21.txt omega -ascii
    elseif (i==22)
        save omegareg22.txt omega -ascii
    elseif (i==23)
        save omegareg23.txt omega -ascii
    elseif (i==24)
        save omegareg24.txt omega -ascii
    elseif (i==25)
        save omegareg25.txt omega -ascii
    elseif (i==26)
        save omegareg26.txt omega -ascii
    elseif (i==27)
        save omegareg27.txt omega -ascii
    elseif (i==28)
        save omegareg28.txt omega -ascii
    elseif (i==29)
        save omegareg29.txt omega -ascii
    elseif (i==30)
        save omegareg30.txt omega -ascii
    elseif (i==31)
        save omegareg31.txt omega -ascii
    elseif (i==32)
        save omegareg32.txt omega -ascii
    elseif (i==33)
        save omegareg33.txt omega -ascii
    elseif (i==34)
        save omegareg34.txt omega -ascii
    elseif (i==35)
        save omegareg35.txt omega -ascii
    elseif (i==36)
        save omegareg36.txt omega -ascii
    elseif (i==37)
        save omegareg37.txt omega -ascii
    elseif (i==38)
        save omegareg38.txt omega -ascii
    elseif (i==39)
        save omegareg39.txt omega -ascii
    elseif (i==40)
        save omegareg40.txt omega -ascii
    elseif (i==41)
        save omegareg41.txt omega -ascii
    elseif (i==42)
        save omegareg42.txt omega -ascii
    elseif (i==43)
        save omegareg43.txt omega -ascii
    elseif (i==44)
        save omegareg44.txt omega -ascii
    elseif (i==45)
        save omegareg45.txt omega -ascii
    elseif (i==46)
        save omegareg46.txt omega -ascii
    elseif (i==47)
        save omegareg47.txt omega -ascii
    end
end
size(ret)
save aggreturns.txt ret -ascii
save aggreturnsp.txt retp -ascii
