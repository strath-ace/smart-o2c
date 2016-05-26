close all
clear all
clc

options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-6);

b0_old =   89.783015631382568*pi/180;
bf_old = -89.620842445292510*pi/180;
Tf_old =   120;%1.090798003571203e+02;
a = 4e-3;
g = 1.6e-3;

sol = fmincon(@(x) x(3),[b0_old bf_old Tf_old],[],[],[],[],[-pi/2 -pi/2 120],[pi/2 pi/2 300],@ascent_constr,options);

b0 = sol(1);
bf = sol(2);
Tf = sol(3);

c = (tan(b0)-tan(bf))/Tf;

t = linspace(0,Tf,1000);

b = atan(tan(b0)-c*t);
u = a/c*log((tan(b0)+sec(b0))./(tan(b)+sec(b)));
v = a/c*(sec(b0)-sec(b))-g*t;
x = a/c^2*( sec(b0)-sec(b)-tan(b).*log( (tan(b0)+sec(b0))./(tan(b)+sec(b)) ) );
y = a/(2*c^2)*( (tan(b0)-tan(b))*sec(b0) - (sec(b0)-sec(b)).*tan(b) - log( (tan(b0)+sec(b0))./(tan(b)+sec(b))) )-0.5*g*t.^2;

figure(1)
plot(t,x)
axis([0 250 0 120])
xlabel('t')
ylabel('x')
figure(2)
plot(t,y)
axis([0 250 0 10])
xlabel('t')
ylabel('y')
figure(3)
plot(t,u)
axis([0 250 0 1])
xlabel('t')
ylabel('u')
figure(4)
plot(t,v)
axis([0 250 0 0.25])
xlabel('t')
ylabel('v')
figure(5)
plot(t,b)
axis([0 250 -pi/2 pi/2])
xlabel('t')
ylabel('\beta')
