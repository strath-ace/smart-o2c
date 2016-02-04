function [f,c]=ZDtestfun(x,mof)
%
%  [f,c]=ZDtestfun(x,mof);
% 
%  Zitzler-Thiele-Deb 2000 test functions for multiobjective optimisation
%
%  INPUT
%           x : state vector
%           mof : funtion name
%
%  OUTPUT 
%           f : multiobjective vector
%

% (c) Massimiliano Vasile  2004
% Modified by Lorenzo A. Ricciardi 2015 (added Fonseca et Al. Function)
c=0;

switch mof
    case 'deb'
        alpha=2;
        q=4;
        f(1)=x(1);
        f(2)=(1+10*x(2))*(1-(x(1)/(1+10*x(2)))^alpha-x(1)*sin(2*pi*q*x(1))/(1+10*x(2)) );
    case 'deb2'
        if (x(1)<0)
            x(1)=-x(1);
        end
        f(1)=x(1);
        g=11+x(2)^2-10*cos(2*pi*x(2));
        
        if (f<g)
            h=1-sqrt(f(1)/g);
        else
            h=0;
        end    
        f(2)=g*h;
    case 'scha'
        if (x(1)<=1)
            f(1)=-x(1);
        elseif(x(1)>1)&(x(1)<=3)
            f(1)=-2+x(1);
        elseif(x(1)>3)&(x(1)<=4)
            f(1)=4-x(1);
        elseif(x(1)>4)
            f(1)=-4+x(1);
        end
        f(2)=(x(1)-5)^2;
     
    case 'kur1'
        n=length(x);
        f(1)=0;
        f(2)=0;        
        for i=1:n-1
            f(1)=f(1)-10*exp(-0.2*sqrt(x(i)^2+x(i+1)^2));
            f(2)=f(2)+abs(x(i))^0.8+5*sin(x(i)^3);
        end
        f(2)=f(2)+abs(x(n))^0.8+5*sin(x(n)^3);
    case 'kur2'
        n=length(x);
        f(1)=0;
        f(2)=0;        
        for i=1:n-1
            f(1)=f(1)+(abs(x(i))^0.8+5*sin(x(i))^3+3.5828);
            f(2)=f(2)+(1-exp(-0.2*sqrt(x(i)^2+x(i+1)^2)));
        end
        f(1)=f(1)+(abs(x(n))^0.8+5*sin(x(n))^3+3.5828);
        
    case 'zdt1'
        n=length(x);
        
        f(1)=abs(x(1));
        g=0;
        for i=2:n
            g=g+abs(x(i))/(n-1);
        end
        g=1+9*g;
        f(2)=g*(1-sqrt(f(1)/g));
    case 'zdt2'
        n=length(x);
        
        f(1)=x(1);
        g=1;
        g2=0;
        for i=2:n
            g2=g2+abs(x(i));
        end
        g=g+9*(g2/(n-1));
        f(2)=g*(1-(f(1)/g)^2);
        
    case 'zdt3'
        n=length(x);
        f(1)=abs(x(1));
        g=1;
        for i=2:n
            g=g+9*abs(x(i))/(n-1);
        end
        
        f(2)=g*(1-sqrt(f(1)/g)-(f(1)/g)*sin(10*pi*f(1)));
        
    case 'zdt4'
        n=length(x);
        f(1)=abs(x(1));
        g=1+10*(n-1);
        for i=2:n
            g=g+x(i)^2-10*cos(4*pi*x(i));
        end
        
        f(2)=g*(1-sqrt(abs(x(1))/g));
        
    case 'zdt6'
        n=length(x);
        f(1)=1-exp(-4*x(1))*sin(4*pi*x(1))^6;
        g=1;
        g2=0;
        for i=2:n
            g2=g2+abs(x(i));
        end
        g=g+9*(g2/(n-1))^0.25;
        
        f(2)=g*(1-(f(1)/g)^2);
        
    case 'quad'
        f(1)=x(1)^2+x(2)^2;
        f(2)=(x(1)-0.5)^2+(x(2)-0.5)^2;    
    case 'constr'
        n=length(x);
        
        f(1)=abs(x(1));
        g=0;
        for i=2:n
            g=g+abs(x(i))/(n-1);
        end
        g=1+9*g;
        f(2)=g*(1-sqrt(f(1)/g));
        the=-0.2*pi;
        a=0.2;
        b=10;
        c=1;
        d=6;
        e=1;
        
        c=cos(the)*(f(2)-e)-sin(the)*f(1)-a*abs(  sin(  b*pi*(  sin(the)*(f(2)-e)+cos(the)*f(1)  )^c) )^d;
        c=-c;
    case 'dtlz6'
        k=20;
        M=3;
        n=M+k-1;
        for i=1:M-1
            f(i)=x(i);
        end
        g=1;
        for i=3:2+k
            g=g+9*x(i)/k;
        end
        h=M;
        for i=1:M-1
            h=h-f(i)*(1+sin(3*pi*f(i)))/(1+g);
        end
        f(M)=(1+g)*h;
    case 'sch'
        f(1) = x^2;
        f(2) = (x-2)^2;        
    case 'fon'
        if length(x)~=3
            error('fon fun must have exactly 3 parameters')
        else
            f(1) = 1-exp(-(sum((x-1/sqrt(3)).^2)));
            f(2) = 1-exp(-(sum((x+1/sqrt(3)).^2)));
        end
    case 'pol'
        A = @(x1,x2) 0.5*sin(x1)-2*cos(x1)+sin(x2)-1.5*cos(x2);
        B = @(x1,x2) 1.5*sin(x1)-cos(x1)+2*sin(x2)-0.5*cos(x2);
        f(1) = 1 + ( A(1,2) - A(x(1),x(2)) )^2 + ( B(1,2) - B(x(1),x(2)));
        f(2) = (x(1)+3)^2+(x(2)+1)^2;
    case 'dtlz1'
        g = 100*(size(x,1)-2)+100*sum( (x(3:end)-0.5).^2 -cos(20*pi*(x(3:end)-0.5)) );
        f(1) = (1+g)*x(1)*x(2);
        f(2) = (1+g)*x(1)*(1-x(2));
        f(3) = (1+g)*(1-x(1));   
    case 'dtlz2'
        g = sum(x(3:end).^2);
        f(1) = (1+g)*cos(x(1)*pi/2)*cos(x(2)*pi/2);
        f(2) = (1+g)*cos(x(1)*pi/2)*sin(x(2)*pi/2);
        f(3) = (1+g)*sin(x(1)*pi/2);   
    otherwise        
        return
end


return

