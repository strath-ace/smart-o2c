function [x, fval, exitflag, output] = optimise_mpaidea(fitnessfcn, LB, UB, options)

%% optimise_mpaidea
% MP-AIDEA-ALR: Multi-Population Adaptive Inflationary Differential Evolution Algorithm with Adaptive Local Restart.
% - Adaptation of CR and F inside the single populations
% - Adaptation of the dimension of the bubble for the local restart
% - Adaptive local/global restart
%
%% Inputs:
%
% * fitnessfcn : function handle to cost function (real function)
% * LB : decision variables lower boundaries
% * UB : decision variables upper boundaries
% * options : structure containing MP-AIDEA specific input information (if empty default values are used), namely
%           * options.population: initial population matrix pop is a 3D matrix where the third dimension identifies the numbers of populations
%           * options.delta_local: dimension of the bubble for the local restart of population. If empty, the value 
%                                  is adapted by MP-AIDEA. If a value is defined for options.delta_local, the
%                                  adaptive mechanism of delta local of MP-AIDEA is not activated
%           * options.delta_global: characteristic dimension for the global restart of the population. delta_global is the distance from the clusters of
%                                   the local minima at which the new elements of the population need to be.
%                                   This parameter is not adapted, therefore it must be defind by the user
%           * options.rho: contraction threshold for the population, that identifies when the Differential Evolution has to be stopped. 
%                          This parameter is not adapted, therefore it must be defined by the user
%           * options.DE_strategy: 1/2, to select the DE strategy. 
%                                  1: DE/Rand and DE/CurrentToBest
%                                  2: DE/Rand and DE/Best
%           * options.prob_DE_strategy: Probability of having one of the two DE strategies defined in options.DE_strategy. If options.DE_strategy is 2,
%                                       options.prob_DE_strategy is the probability of using DE/Rand rather than DE/Best
%           * options.dd_CRF: parameter for the adaptation of CRF
%           * options.nFeValMax: maximum number of function evaluations
%
%% Output:
% * x : optimal design solution found
% * fval: optimal cost
% * exitflag: exitflag from fmincon
% * output: structure containing MP-IDEA specific ouput information, namely
%           * output.memories: archive containing all the local minima plus the population for each restart. The solutions are sorted from the best to the worst. One cell for each population
%           * output.archivebest: archive of local minima (obtained by local search)
%           * output.population_evolution: initial and final population obtained at each DE step for each population
%           * output.vval_evolution: for each population, important DE parameters 
%           * output.B_mean: mean value of the matrix for the adaptation of delta local
%           * output.delta_local: value of delta_local for each population
%           * output.number_LR: number of local restart for each population
%           * output.number_GR: number of global restart for each population
%
%% Author(s): Edmondo Minisci and Massimiliano Vasile (2013), Marilena Di Carlo (2015)
% email: edmondo.minisci@strath.ac.uk massimiliano.vasile@strath.ac.uk marilena.di-carlo@strath.ac.uk
%
%% References:
% * M. Di Carlo, M. Vasile, E. Minisci, "Multi-Population Adaptive Inflationary Differential Evolution Algorithm with Adaptive Local Restart", IEEE Congress on Evolutionary Computation, CEC 2015, Sendai, Japan
% https://pure.strath.ac.uk/portal/en/publications/multipopulation-adaptive-inflationary-differential-evolution-algorithm-with-adaptive-local-restart(a3a478f0-54ad-4bc0-95d8-30bc45fb3e2f).html

% * M. Di Carlo, M. Vasile, E. Minisci, "Multi-Population Adaptive Inflationary Differential Evolution Algorithm", 2014 BIOMA (Bio-Inspired % Optimisation Methods and their Applications) Workshop, Ljubljana, Slovenia
% https://pure.strath.ac.uk/portal/en/publications/multipopulation-adapative-inflationary-differential-evolution(fef0a381-e18e-43f4-ae50-db261f0f755d).html

% * E. Minisci, M. Vasile "Adaptive Inflationary Differential Evolution Algorithm", IEEE Congress on Evolutionary Computation, CEC 2014, China
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6900587


% Run MP-AIDEA
[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, LB, UB, options.population, options);

% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), numel(LB));
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = memories_record(end, end, i);
end

output.memories             = memories;
output.archivebest          = archivebest;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;

end
