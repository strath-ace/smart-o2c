close all
clear all
clc

%% load files

load('Ascent/Ascent_multiple_runs_10k_each_standard_pattern_equidistant_reposit_6th_order_10pts/10_run_memory_1.mat')

mm = [];

for i=1:10
    
    mm = [mm;mem(i).memory];
    
end

load('Ascent/Ascent_multiple_runs_10k_each_standard_pattern_equidistant_reposit_6th_order_10pts/10_run_memory_2.mat')

for i=1:10
    
    mm = [mm;mem(i).memory];
    
end

load('Ascent/Ascent_multiple_runs_10k_each_standard_pattern_equidistant_reposit_6th_order_10pts/10_run_memory_3.mat')

for i=1:10
    
    mm = [mm;mem(i).memory];
    
end

clear mem

%% extract non dominated front

dom=dominance(mm(:,end-3:end-2),0);
nondom = mm(dom==0,:);

%% load gradient solutions

load('C:\Users\LA\Dropbox\SMART\smart-o2c\Problems\MACS\MO_control\Ascent\Ascent_gradient_500kfuneval\10_run_memory_1.mat')

mm2 = [];

for i=1:10
    
    mm2 = [mm2;mem(i).memory];
    
end

load('C:\Users\LA\Dropbox\SMART\smart-o2c\Problems\MACS\MO_control\Ascent\Ascent_gradient_500kfuneval\10_run_memory_2.mat')

for i=1:10
    
    mm2 = [mm2;mem(i).memory];
    
end

load('C:\Users\LA\Dropbox\SMART\smart-o2c\Problems\MACS\MO_control\Ascent\Ascent_gradient_500kfuneval\10_run_memory_3.mat')

for i=1:10
    
    mm2 = [mm2;mem(i).memory];
    
end

clear mem

%% extract non dominated front of gradient solutions

dom=dominance(mm2(:,end-3:end-2),0);
nondom2 = mm2(dom==0,:);


%% select 4 uniformly distributed points

mins = min(nondom(:,end-3:end-2));
maxs = max(nondom(:,end-3:end-2));

[selection,dd,energy,ener2,mins,maxs]=arch_shrk6([],[],nondom(1:4,:),0,[],mins,maxs,size(nondom,2)-4,2,4);
[selection,dd,energy,ener2,mins,maxs]=arch_shrk6(selection,dd,nondom(5:end,:),energy,ener2,mins,maxs,size(nondom,2)-4,2,4);

%% sort chosen points

[~,b] = sort(selection(:,end-3));

qq = selection(b,:);    %sort wrt t_f

%% select 4 points from gradient solutions, with closest mission time to the previous ones

qq2 = zeros(size(qq));

for i = 1:size(qq,1)
   
    [dist,index] = min(abs(nondom2(:,end-3)-qq(i,end-3)));
    qq2(i,:) = nondom2(index,:);
    
end

%% Plot Pareto front

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 10)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 10)

plot(nondom(:,end-3),nondom(:,end-2),'b.')
box on
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(nondom2(:,end-3),nondom2(:,end-2),'g.')
plot(qq(:,end-3),qq(:,end-2),'r+','MarkerSize', 10)
plot(qq2(:,end-3),qq2(:,end-2),'m+','MarkerSize', 10)
xlabel('t_f')
ylabel('-v_x(t_f)')
% text(qq(1,end-3)+2,qq(1,end-2)-0.02,'1')
% text(qq(2,end-3)+2,qq(2,end-2)+0.02,'2')
% text(qq(3,end-3)+2,qq(3,end-2)+0.02,'3')
% text(qq(4,end-3),qq(4,end-2)+0.02,'4')
%text(qq(5,end-3),qq(5,end-2)+0.02,'5')

%% Compute metrics wrt Pareto front

nonnorm = (nondom(:,end-3:end-2)-repmat(min(nondom(:,end-3:end-2)),size(nondom,1),1));%./repmat(max(nondom(:,end-3:end-2)),size(nondom,1),1);
nonnorm = nonnorm./repmat(max(nonnorm(:,1:2)),size(nonnorm,1),1);

%fp = nondom(:,end-3:end-2);
fp = nonnorm;
xp = nondom(:,1:end-4);

M1 = zeros(30,1);
M4 = M1;
AveGDP = M1;
AveIGDP = M1;
AveHaus = M1;

weights = [1 1];

% for j=1:30
%
%     this = mm((j-1)*10+1:j*10,:);
%
%     x = this(:,1:end-4);
%     f = this(:,end-3:end-2);
%     f = f-repmat(min(nondom(:,end-3:end-2)),size(f,1),1);
%     f = f./repmat(max(nondom(:,end-3:end-2))-min(nondom(:,end-3:end-2)),size(f,1),1);
%
%     [M1(j),~,~,M4(j)] = ZDTmetrics(x,f,xp,fp,0,weights);
% 	[AveIGDP(j),AveGDP(j),AveHaus(j)]=AveHausMetrics(x,f,xp,fp,inf);
%
% end
%
% M1mean = mean(M1);
% M1var = var(M1);
% M4mean = mean(M4);
% M4var = var(M4);
% AveIGDPmean = mean(AveIGDP);
% AveIGDPvar = var(AveIGDP);
% AveGDPmean = mean(AveGDP);
% AveGDPvar = var(AveGDP);
% AveHausmean = mean(AveHaus);
% AveHausvar = var(AveHaus);

%% DFET PROBLEM TRANSCRIPTION

% Equations, initial conditions and time span
a=4e-3;
f = @(x,u,t) [x(2); a*cos(u(1)); x(4); -0.0016+a*sin(u(1))];
dfx = @(x,u,t) [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
dfu = @(x,u,t) [0; -a*sin(u(1)); 0; a*cos(u(1))];
dft = @(x,u,t) [0; 0; 0; 0];

t_0 = 0;

x_0 = [0 0 0 0];

imposed_final_states = [0 0 1 1];   % mask vector with the variables on which a final condition is imposed
x_f = [0 0.1 10 0];                % vector of final conditions (size???)

% Discretisation settings

num_elems = 4;
state_order = 6;
control_order = 6;
DFET = 1;
state_distrib = 'Lobatto'; % or Cheby
control_distrib = 'Legendre';

integr_order = 2*state_order;%2*state_order-1;

num_eqs = length(x_0);
num_controls = 1;

% Make checks

test_order = state_order+(DFET==1);

if DFET==0
    
    total_num_equations = num_elems*(test_order+1)*num_eqs+sum(imposed_final_states);
    total_num_unknowns = num_elems*(state_order+1)*num_eqs+num_elems*(control_order+1)*num_controls;
    
else
    
    total_num_equations = ((test_order+1)+(test_order)*(num_elems-1))*num_eqs;
    total_num_unknowns = ((test_order)*num_elems)*num_eqs+num_elems*(control_order+1)*num_controls+sum(imposed_final_states==0);
    
end

if total_num_equations>total_num_unknowns
    
    error('Problem is over-constrained: either reduce number of final constraints, increase number of control variables, use higher order polynomials for control variables, or use more elements');
    
end

if total_num_equations==total_num_unknowns
    
    warning('Number of constraints equal to number of unknown control coefficients: optimal control is not possible, only constraints satisfaction');
    
end

% Transcribe constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

state_bounds = [-inf inf;-inf inf; -inf inf;-inf inf];
control_bounds = [-pi/2 pi/2];

[vlb,vub] = transcribe_bounds(state_bounds,control_bounds,structure);

vlb = [100;vlb]';
vub = [250;vub]';

tol_conv = 1e-6;
maxits = 10000;

fminconoptions = optimset('Display','off','MaxFunEvals',maxits,'Tolcon',tol_conv,'GradCon','on');

%% plot trajectories

figure(2)
box on
set(gca,'LooseInset',get(gca,'TightInset'));

for i = 1:size(qq,1)
    
    [x,u,xb] = extract_solution(qq(i,2:length(vlb)),structure,x_f);
    tplot = linspace(t_0,qq(i,1),100);
    [xt,ut] = eval_solution_over_time(x,u,t_0,qq(i,1),tplot,structure.uniform_els,structure);
    plot(xt(:,1),xt(:,3))
    hold on
    
    % gradient solutions
    [x,u,xb] = extract_solution(qq2(i,2:length(vlb)),structure,x_f);
    tplot = linspace(t_0,qq2(i,1),100);
    [xt,ut] = eval_solution_over_time(x,u,t_0,qq2(i,1),tplot,structure.uniform_els,structure);
    plot(xt(:,1),xt(:,3))
    
end
xlabel('x');
ylabel('y');
%legend('Trajectory 1','Trajectory 1(grad)','Trajectory 2','Trajectory 2(grad)','Trajectory 3','Trajectory 3(grad)','Trajectory 4','Trajectory 4(grad)','Location','East')

%% plot velocities and control laws
hs = {};
hc = {};
for i = 1:size(qq,1)
    
    [x,u,xb] = extract_solution(qq(i,2:length(vlb)),structure,x_f);
    tplot = linspace(t_0,qq(i,1),100);
    [xt,ut] = eval_solution_over_time(x,u,t_0,qq(i,1),tplot,structure.uniform_els,structure);
    
    [hs{i},hc{i}] = plot_solution_vs_time2(x,u,x_0,xb,t_0,qq(i,1),structure.uniform_els,structure,2*i+1);
    hold on
    % gradient solution
    
    [x2,u2,xb2] = extract_solution(qq2(i,2:length(vlb)),structure,x_f);
    tplot = linspace(t_0,qq2(i,1),100);
    [xt2,ut2] = eval_solution_over_time(x2,u2,t_0,qq2(i,1),tplot,structure.uniform_els,structure);
    
    [hs2{i},hc2{i}] = plot_solution_vs_time2(x2,u2,x_0,xb2,t_0,qq2(i,1),structure.uniform_els,structure,2*i+1);
    figure(2*i+1)
    
    % compute analytic solution
    
    options = optimset('Display','off','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-6);
    
    b0_old =   89.783015631382568*pi/180;
    bf_old = -89.620842445292510*pi/180;
    Tf_old =   qq(i,1);
    a = 4e-3;
    g = 1.6e-3;
    
    sol = fmincon(@(x) x(3),[b0_old bf_old Tf_old],[],[],[],[],[-pi/2 -pi/2 qq(i,1)],[pi/2 pi/2 300],@ascent_constr,options);
    
    b0 = sol(1);
    bf = sol(2);
    Tf = sol(3);
    
    c = (tan(b0)-tan(bf))/Tf;
    
    t = linspace(0,Tf,100);
    
    b = atan(tan(b0)-c*t);
    u_analytic = a/c*log((tan(b0)+sec(b0))./(tan(b)+sec(b)));
    v_analytic = a/c*(sec(b0)-sec(b))-g*t;
    x_analytic = a/c^2*( sec(b0)-sec(b)-tan(b).*log( (tan(b0)+sec(b0))./(tan(b)+sec(b)) ) );
    y_analytic = a/(2*c^2)*( (tan(b0)-tan(b))*sec(b0) - (sec(b0)-sec(b)).*tan(b) - log( (tan(b0)+sec(b0))./(tan(b)+sec(b))) )-0.5*g*t.^2;
    
    clf
    plot(tplot,xt(:,2),'b')
    box on
    set(gca,'LooseInset',get(gca,'TightInset'));
    hold on
    plot(tplot,xt2(:,2),'g')
    plot(tplot,xt2(:,4),'b--')
    plot(tplot,xt2(:,4),'g--')    
    plot(t,u_analytic,'b.-');
    plot(t,v_analytic,'b:');
    axis([0 250 0 1])
    xlabel('t')
    ylabel('velocities')
    legend('v_x','v_x(grad)','v_y','v_y(grad)','v_x(analytic)','v_y(analytic)')
    figure(2*i+2)
    box on
    set(gca,'LooseInset',get(gca,'TightInset'));
    hold on
    hnew = plot(t,b,'b--');
    axis([0 250 -pi pi])
    xlabel('t')
    ylabel('controls')
    legend([hc{i}(1),hc2{i}(1),hnew],'u(MACS)','u(gradient)','u(analytic)')
    figure(1)
    plot(t(end),-u_analytic(end),'xk')
    drawnow
    
end 

%% compare extreme points with single objective solutions

% compute analytic solution

options = optimset('Display','off','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-6);

b0_old =   89.783015631382568*pi/180;
bf_old = -89.620842445292510*pi/180;
Tf_old =   qq(1,1);
a = 4e-3;
g = 1.6e-3;

sol = fmincon(@(x) x(3),[b0_old bf_old Tf_old],[],[],[],[],[-pi/2 -pi/2 Tf_old],[pi/2 pi/2 300],@ascent_constr,options);

b0 = sol(1);
bf = sol(2);
Tf = sol(3);

c = (tan(b0)-tan(bf))/Tf;

t = linspace(0,Tf,100);

b = atan(tan(b0)-c*t);
u_analytic = a/c*log((tan(b0)+sec(b0))./(tan(b)+sec(b)));
v_analytic = a/c*(sec(b0)-sec(b))-g*t;
x_analytic = a/c^2*( sec(b0)-sec(b)-tan(b).*log( (tan(b0)+sec(b0))./(tan(b)+sec(b)) ) );
y_analytic = a/(2*c^2)*( (tan(b0)-tan(b))*sec(b0) - (sec(b0)-sec(b)).*tan(b) - log( (tan(b0)+sec(b0))./(tan(b)+sec(b))) )-0.5*g*t.^2;

figure(3)
clf
[x,u,xb] = extract_solution(qq(1,2:length(vlb)),structure,x_f);
tplot = linspace(t_0,qq(1,1),10);
[xt,ut] = eval_solution_over_time(x,u,t_0,qq(1,1),tplot,structure.uniform_els,structure);
plot(tplot,xt(:,2),'bo')
box on
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(tplot,xt(:,4),'b*')
axis tight
tplot = linspace(t_0,qq(1,1),100);
load('Ascent/min_time_states.mat')
load('Ascent/min_time_controls.mat')
plot(tplot,xt(:,2),'b')
plot(tplot,xt(:,4),'b--')
plot(t,u_analytic,'b:');
plot(t,v_analytic,'b.-');
xlabel('t')
ylabel('velocities')
legend('v_x(MACS)','v_y(MACS)','v_x(single obj)','v_y(single obj)','v_x(analytic)','v_y(analytic)')
axis tight


figure(4)
clf
[x,u,xb] = extract_solution(qq(1,2:length(vlb)),structure,x_f);
tplot = linspace(t_0,qq(1,1),10);
[xt,ut] = eval_solution_over_time(x,u,t_0,qq(1,1),tplot,structure.uniform_els,structure);
plot(tplot,ut,'b*')
box on
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
axis tight
tplot = linspace(t_0,qq(1,1),100);
load('Ascent/min_time_states.mat')
load('Ascent/min_time_controls.mat')
plot(tplot,ut,'b')
plot(t,b,'b:')
xlabel('t')
ylabel('controls')
legend('u(MACS)','u(single obj)','u(analytic)')
axis tight

% compute analytic solution

options = optimset('Display','off','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-6);

b0_old =   89.783015631382568*pi/180;
bf_old = -89.620842445292510*pi/180;
Tf_old =   qq(end,1);
a = 4e-3;
g = 1.6e-3;

sol = fmincon(@(x) x(3),[b0_old bf_old Tf_old],[],[],[],[],[-pi/2 -pi/2 qq(end,1)],[pi/2 pi/2 300],@ascent_constr,options);

b0 = sol(1);
bf = sol(2);
Tf = sol(3);

c = (tan(b0)-tan(bf))/Tf;

t = linspace(0,Tf,100);

b = atan(tan(b0)-c*t);
u_analytic = a/c*log((tan(b0)+sec(b0))./(tan(b)+sec(b)));
v_analytic = a/c*(sec(b0)-sec(b))-g*t;
x_analytic = a/c^2*( sec(b0)-sec(b)-tan(b).*log( (tan(b0)+sec(b0))./(tan(b)+sec(b)) ) );
y_analytic = a/(2*c^2)*( (tan(b0)-tan(b))*sec(b0) - (sec(b0)-sec(b)).*tan(b) - log( (tan(b0)+sec(b0))./(tan(b)+sec(b))) )-0.5*g*t.^2;

figure(9)
clf
[x,u,xb] = extract_solution(qq(end,2:length(vlb)),structure,x_f);
tplot = linspace(t_0,qq(end,1),10);
[xt,ut] = eval_solution_over_time(x,u,t_0,qq(end,1),tplot,structure.uniform_els,structure);
plot(tplot,xt(:,2),'bo')
box on
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(tplot,xt(:,4),'b*')
axis ([0 250 0 1])
tplot = linspace(t_0,qq(end,1),100);
load('Ascent/max_vel_states.mat')
load('Ascent/max_vel_controls.mat')
plot(tplot,xt(:,2),'b')
plot(tplot,xt(:,4),'b--')
plot(t,u_analytic,'b:');
plot(t,v_analytic,'b.-');
xlabel('t')
ylabel('velocities')
legend('v_x(MACS)','v_y(MACS)','v_x(single obj)','v_y(single obj)','v_x(analytic)','v_y(analytic)')

figure(10)
% clf
% [x,u,xb] = extract_solution(qq(end,2:length(vlb)),structure,x_f);
% tplot = linspace(t_0,qq(end,1),10);
% [xt,ut] = eval_solution_over_time(x,u,t_0,qq(end,1),tplot,structure);
% plot(tplot,ut,'b*')
% box on
% set(gca,'LooseInset',get(gca,'TightInset'));
% hold on
% axis ([0 250 -pi pi])
tplot = linspace(t_0,qq(end,1),100);
load('Ascent/max_vel_states.mat')
load('Ascent/max_vel_controls.mat')
hnew = plot(tplot,ut,'b--');
hnew2 = plot(t,b,'b:');
xlabel('t')
ylabel('controls')
legend([hc{end}(1),hnew,hnew2],'u(MACS)','u(single obj)','u(analytic)')
