function [M1,M2,M3,M4]=ZDTmetrics(x,f,xp,fp,sigma,weights);
%
%   [M1,M2,M3,M4]=ZDTmetrics(x,f,xp,fp,sigma,weights)
%
%   Zitzler & Thiele Convergence Metrics
%
%   INPUT
%               x     : achieved solution set
%               f     : achieved criteria set
%               xp    : true Pareto solution set
%               fp    : true Pareto front
%               sigma : neighborhood parameter
%               weihts: weights on f components
%
%   OUTPUT
%               M1 : front convergence metrics
%               M2 : distribution metrics
%               M3 : Front metrics
%               M4 : local convergence
%
%  from " Comparison of Multiobjective Evolutionary Algorithms: Empirical
%  Results", Zitzler, Thiele, Deb, Evolutionary Computation 2000

% (c) Massimiliano Vasile 2005

lx=length(x(:,1));
lf=length(f(:,1));
lxp=length(xp(:,1));
lfp=length(fp(:,1));
for i=1:lfp
fp(i,:)=fp(i,:).*weights;
end
for i=1:lf
f(i,:)=f(i,:).*weights;
end

% Spreading metrics
M1=0;
for i=1:lfp
    dfp=1e36;
    for j=1:lf
        if all(isnan(f(j,:)))<1            
%            dfp=min([dfp,100*abs(norm((f(j,:)-fp(i,:))./(1e-8+fp(i,:))))]);
            dfp=min([dfp,norm(f(j,:)-fp(i,:))]);
            
        end
    end
   
    M1=M1+dfp;
end

M1=M1/lfp;

% Convergence metrics
M4=0;
for i=1:lf
    dfp=1e36;
    for j=1:lfp
        if all(isnan(f(i,:)))<1
%            dfp=min([dfp,abs(norm((f(i,:)-fp(j,:))./(1e-8+fp(i,:))))]);
            dfp=min([dfp,abs(norm((f(i,:)-fp(j,:))))]);

        end
    end
    M4=M4+dfp;
end

M4=M4/lf;

% Distribution metrics
M2=0;
for i=1:lf
    for j=1:lf
        if norm(f(i,:)-f(j,:))>sigma
            M2=M2+1;
        end
    end
end
M2=M2/(lf-1);

% Front metrics
M3=0;
for i=1:lf
    dfp=0;
    for j=1:lf
        dfp=max([dfp,norm(f(i,:)-f(j,:))]);
    end
    M3=M3+dfp;
end
M3=sqrt(M3);

return