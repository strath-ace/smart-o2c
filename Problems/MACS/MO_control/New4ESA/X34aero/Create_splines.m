close all
clear
clc

load('alpha');
load('Mach');
load('Cl');
load('Cd');

% autobuild function :D

load('cl_model_values');
load('cl_model_names');

for i = 1:length(cl_model_values)

    eval([cl_model_names{i},'=',num2str(cl_model_values(i)),';']);
    
end

cl_fun = @(x,y) (a0+a1*x+a2*x^2)+(b0+b1*x+b2*x^2)*(k1/l1)*((y-s1)/l1)^(k1-1)*exp(-((y-s1)/l1)^k1)+(c0+c1*x+c2*x^2)*(k2/l2)*((y-s3)/l2)^(k2-1)*exp(-((y-s3)/l2)^k2);

% autobuild function :D

load('cd_model_values');
load('cd_model_names');

for i = 1:length(cd_model_values)

    eval([cd_model_names{i},'=',num2str(cd_model_values(i)),';']);
    
end

cd_fun = @(x,y) (a0+a1*x+a2*x^2+a3*x^3)+(b0+b1*x+b2*x^2+b3*x^3)*(k1/l1)*((y-s1)/l1)^(k1-1)*exp(-((y-s1)/l1)^k1)+(c0+c1*x+c2*x^2+c3*x^3)*(k2/l2)*((y-s3)/l2)^(k2-1)*exp(-((y-s3)/l2)^k2);

[AoA,M] = meshgrid(alpha,Mach);

aa = AoA(:)';
mm = M(:)';
ll = Cl(:)';
dd = Cd(:)';

alpha_plot = 0:0.1:35;
Mach_plot = 0:0.1:9;

[Ap,Mp] = meshgrid(alpha_plot,Mach_plot);

Clp = 0*Ap;
Cdp = 0*Ap;

for i = 1:size(Ap,1)
   
    for j = 1:size(Ap,2)
       
        Clp(i,j) = cl_fun(Ap(i,j),Mp(i,j));
        Cdp(i,j) = cd_fun(Ap(i,j),Mp(i,j));

    end
    
end

figure(1)
ss = surf(Mp,Ap,Clp);
ss.LineStyle = 'None';
hold on
plot3(M,AoA,Cl,'k*')
xlabel('\alpha')
ylabel('Mach')
zlabel('Cl')

figure(2)
ss2 = surf(Mp,Ap,Cdp);
ss2.LineStyle = 'None';
hold on
plot3(M,AoA,Cd,'k*')
xlabel('\alpha')
ylabel('Mach')
zlabel('Cd')

figure(3)
ss3 = surf(Mp,Ap,Clp./Cdp);
ss3.LineStyle = 'None';
hold on
plot3(M,AoA,Cl./Cd,'k*')
xlabel('\alpha')
ylabel('Mach')
zlabel('E')