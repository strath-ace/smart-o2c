function [AveIGD,AveGD,AveHaus]=AveHausMetrics(x,f,xp,fp,exp)
%
%   [AveIGD,AveGD,AveHaus]=AveHausMetrics(x,f,xp,fp)
%
%   Oliver Schutze, Xavier Esquivel, Adriana Lara, and Carlos A. Coello Coello
%
%   INPUT
%               x     : achieved solution set
%               f     : achieved criteria set
%               xp    : true Pareto solution set
%               fp    : true Pareto front
%               exp   : exponent used for norm
%
%   OUTPUT
%               AveIGD : corrected IGD
%               AveGD  : corrected GD
%               AveHaus : averaged Hausdorff metric
%
%  from "Measuring the Averaged Hausdorff Distance to the Pareto Front of a
%  Multi-Objective Optimization Problem
%
% (c) Lorenzo Ricciardi 2015

lx=length(x(:,1));
lf=length(f(:,1));
lxp=length(xp(:,1));
lfp=length(fp(:,1));

% Corrected IGD

if exp<inf

    AveIGD=0;
    
else
   
    AveIGD=-1;
    
end

for i=1:lfp
    
    dfp=inf;
    
    for j=1:lf
        
        if all(isnan(f(j,:)))<1            

            dfp=min([dfp,norm(f(j,:)-fp(i,:))]);
            
        end
        
    end
    
    if exp<inf
   
        AveIGD=AveIGD+dfp^exp;
        
    else
       
        AveIGD=max([AveIGD,dfp]);
        
    end
    
end

if exp<inf
    
    AveIGD=(AveIGD/lfp)^(1/exp);

end

% Corrected GD

if exp<inf

    AveGD=0;

else
    
    AveGD=-1;
    
end

for i=1:lf
    
    dfp=inf;
    
    for j=1:lfp
        
        if all(isnan(f(i,:)))<1

            dfp=min([dfp,norm(f(i,:)-fp(j,:))]);

        end
        
    end
    
    if exp<inf
    
        AveGD=AveGD+dfp^exp;
    
    else
       
        AveGD=max([AveGD,dfp]);

    end
    
end

if exp<inf

    AveGD=(AveGD/lf)^(1/exp);

end

AveHaus = max([AveIGD,AveGD]);

return