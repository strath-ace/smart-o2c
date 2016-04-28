function [f,f2,f3,f4,f5,f6,f7]=moo_robust_test(y,par,flag)

nx=length(y)/2;
u=y(nx+1:end);
x=y(1:nx);
mootestfun='simone1';
switch mootestfun
    case 'MVMOO1'
        b=linspace(0.5,nx-1,nx);
        f(1)=sum(x.^2)/nx+0.1*(prod(abs(-u-b)))^(2/nx);
        f(2)=sqrt(30-f(1))+0.1*((u(1)-nx)-sum(u(2:end)-linspace(nx-1,1,nx-1)))^2/nx^2;
    case 'MVMOO2'
        b=linspace(0.5,nx-1,nx);
        f(1)=sum(x.^2)/nx+0.1*(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=sqrt(30-f(1))+0.1*((u(1)-nx)*cos(2*u(1))-sum((u(2:end)-linspace(nx-1,1,nx-1)).*cos(3*u(2:end))))^2/nx^2;
    case 'MVMOO3'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(0.5,nx-1,nx);
        f(1)=(sum(x.^2)/nx+1)*(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=sqrt(30-sum(x.^2)/nx)-fA;
    case 'MVMOO4'

        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(0.5,nx-1,nx);
        
        f(1)=(sum(x.^2)/nx+1)*(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=sqrt(30-sum(x.^2)/nx)*fA;

    case 'MVMOO5'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(0.5,nx-1,nx);
        f(1)=(sum(x.^2)/nx+1)*(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=1/(1+f(1))-fA;
    case 'MVMOO6'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(0.5,nx-1,nx);

        s = sum(-x.*sin(5*2*pi*sqrt(abs(x))));
        y = 5+s/nu;
        
        f(1)=y+(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=1/(1+f(1))-fA;
    case 'MVMOO7'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(1,nx-1,nx);

        s = sum(-x.*sin(5*2*pi*sqrt(abs(x))));
        y = 5+s/nu;
        
        f(1)=y*(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=sqrt(50-y)-fA;
 case 'MVMOO8'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(1,nx-1,nx);

        s = sum(-x.*sin(5*2*pi*sqrt(abs(x))));
        y = 5+s/nu;
        
        f(1)=y+(prod(abs((-u-b).*cos(2*u))))^(2/nx);
        f(2)=sqrt(50-f(1))-fA;
case 'MVMOO9'
        sum1=0;
        sum2=0;
       
        nu=length(u);
        e=4.481689070338064e+000;
        for i=1:nu
            sum1=sum1+(u(i))^2;
            sum2=sum2+1.5*cos(2*pi*2*(u(i)));
        end
        
        fA=-20*exp(-0.2*sqrt(sum1/nu))-exp(sum2/nu)+20+e;
        
        b=linspace(1,nx-1,nx);

        s = sum(-x.*sin(5*2*pi*sqrt(abs(x))));
        y = 5+s/nu;
        
        f(1)=y+(prod(abs((-u-x).*cos(2*u))))^(2/nx);
        f(2)=sqrt(50-f(1))-fA;

        
    case 'simone1'
        % Multi objective test function in d: ZDT2 (known Pareto front)
 
        n = numel(x);
        g1 = x(1);
        q = 1 + 9/(n-1)*sum(x(2:n));
        r = 1 - (g1/q)^2;
        g2 = q*r;
        
        % Single objective test functions in u
        m = numel(u);
        h1 = ( 10*m + sum(u.^2 - 10*cos(2*pi*u)) );                 % Rastrigin function, f(0,...,0) = 0
        h2 = ( 418.9829*m - sum( 100*u.*sin(sqrt(abs(100*u))) ));   % Schwefel function, f(4.209687,...,4.209687) = 0
      
        % Resulting test function
        f(1) = g1 + h1;
        f(2) = g2 + h2;
        
end
f2=0;
f3=0;
f4=0;
f5=0;
f6=0;
f7=0;
return