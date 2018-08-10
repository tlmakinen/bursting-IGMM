function error = errorcorrectedtwostate4(p,d3,pon,dt,doplot,N,L,lesvar,tPol,color)

if ~exist('tPol','var')
    tPol=6;
end
if numel(p)==1
    p=[p 1];    % No normalization
end

LL3=L;
aa=1-exp(p(1));
timelength=size(d3,1);
timepoints=[0:(timelength-1)]*dt;
kon=pon*aa;
koff=aa-kon;

 %%  Generate the L(i), L(j)  3prime
U=[-koff kon; koff -kon  ];
lightperiod=max(size(L));
X0=[1 0]';

[V,D]=eig(U,'vector');
invV=inv(V);
cor3=zeros(1,N);
% Create a look up table for expm(U*((s-1)*dt/tPol-i+j))
expU={};
% Calculate the last of the unconnected autocorrelation function
for s=N
    [expU,order]=expm_4(V,D,invV,(s-1)*dt/tPol,lightperiod-1);
    for i=1:lightperiod
        for j=1:lightperiod
            cor3(1,s)=cor3(1,s)+LL3(i)*LL3(j)*(pon*X0'*expU{order(i-j+lightperiod)}*X0);
        end
    end
    cor3(1,s)=real(cor3(1,s));
end
% Calculate the remaining point of the autocorrelation function
    skip=0;
for s=1:N-1
    if s>2
        if (abs(cor3(1,s-1)-cor3(1,s-2))/cor3(1,N)<1e-4)&&(abs(cor3(1,s-1)-cor3(1,N))/cor3(1,N)<1e-4)
            skip=1;
        end
    end
    if skip
        cor3(1,s)=cor3(1,N);
    else
        [expU,order]=expm_4(V,D,invV,(s-1)*dt/tPol,lightperiod-1);
        for i=1:lightperiod
            for j=1:lightperiod
                cor3(1,s)=cor3(1,s)+LL3(i)*LL3(j)*(pon*X0'*expU{order(i-j+lightperiod)}*X0);
            end
        end
    end
    cor3(1,s)=real(cor3(1,s));
end


finally=zeros(1,timelength);

A=cor3;

for r=0:(timelength-1)
    finally(r+1)=A(r+1) + (1/N)*(1/N - 2/(N-r))*(N*A(1) + sum(2*(N-[1:(N-1)]).*A(2:N))) +2/(N*(N-r))*(r*A(1) + sum( 2*(r-[1:(r-1)]).*A(2:r)) +sum(A(2:N).*(min(N,r+[1:(N-1)])-max(r,[1:(N-1)]) )      )  ) ;
end


cor3=finally;
% Performing the normalization:
cor3=cor3*p(2)/cor3(2)*d3(2);
%lesvar=1./(timelength:-1:1)';
err3=sum((d3(2:timelength)-cor3(2:timelength)').^2.*(timelength-1:-1:1)');
error=err3;
if isnan(error)
    error=1e20;
end

if doplot>0
    if ~exist('color','var')
        color='r';
    end
    hold on;
    if any(doplot==1)
        plot(timepoints,cor3,color,'Linewidth',1);
    end
    if any(doplot==2)
        shadedErrorBar(timepoints,d3,sqrt(lesvar),'b--',0.5);
    end
    if any(doplot==3)
        plot(timepoints,d3,'b--');
    end
    if any(doplot==4)
        plot(timepoints,d3,'g--');
    end
    hold on;
    plot(timepoints,cor3,color,'Linewidth',1);
    h=xlim;
    h(1)=0;
    xlim(h);
end