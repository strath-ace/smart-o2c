
%mex cec14_func.cpp -DWINDOWS

% Number of populations
pop_number = 4;

% Number of individuals in the population
NP = 10;

% Dimension of the problem
D = 10;

% Lower and upper boundary search space
vub = 100*ones(1,D);
vlb = -100*ones(1,D);

% Options for MP-AIDEA
options.delta_global = 0.1;
options.rho = 0.2;
options.prob_DE_strategy = 0.5;
options.dd_CRF = 3;
options.nFeValMax = 100000;
options.DE_strategy = 1;
options.plot_flag = 0;

% Create the population
population = zeros(NP,D,pop_number);

for s = 1 : pop_number
    pop = lhsdesign(NP,D,'criterion','maximin').*repmat(vub-vlb,NP,1)+repmat(vlb,NP,1);
    population(:,:,s) = pop;
end

% Function to optimize
func_num = 1;

input.population = population;
input.vlb = vlb;
input.vub = vub;

 % Optimize   
[memories, B_mean, bubble, archivebest,options] = MP_AIDEA_ALR(@(x)cec14_func(x,func_num),input, options);



