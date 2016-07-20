function [x,fval,exitflag,output] = optimise_macs(fitnessfcn,LB,UB,options,varargin)

%% optimise_macs
%% Inputs:
%
% * fitnessfcn : function handle to cost function (vectorial function)
% * LB : decision variables lower boundaries
% * UB : decision variables upper boundaries
% * options : structure of optimization parameters
%           * options.maxnfeval: max num of function evaluations albeit forced in a "weak" way (default 5000)
%           * options.popsize: max num of agents (default 10)
%           * options.rhoini: initial local hypercube size (default 1)
%           * options.F: F parameter for Differential Evolution (default 0.9)
%           * options.CR: crossover probability (default 0.9)
%           * options.p_social: ratio between population performing social moves and total population (default 0.2)
%           * options.max_arch: max size of the archive at output (default 10)
%           * options.coord_ratio: ratio of coordinates to check after which no further search will be performed, as in Zuiani 3b, pag 6 (default 0.25)
%           * options.contr_ratio: contraction ratio for rho (default 0.5)
%           * options.draw_flag: plot flag (default 0)
%           * options.cp: constraint flag (default 0)
%           * options.MBHflag: MBH steps flag (default 0) Currently only for unconstrained
%           * options.DE_strategy: strategy to use for DE (ie pull towards best element or random one, default 'best')
% * varargin : additional fitness input parameters (not decision variables)
%
%% Output:
% * x : archive containing the solution vectors, matrix nxd where n is the number of solution in the Pareto frontier and d is the dimension of the decision vector
% * fval : archive containing objective function values, matrix nxm where n is the numbe of solution in the Pareto forntier and m the number of objectives
% * exitflag : always return 1 (unused)
% * output : additional outputs from MACS
%          * output.memory : the entire archive
%          * output.nfeval : effective number of objective function evaluations
%          * output.ener : value of the energy of the archive, used by Energy
%                   Based Archiving strategy
%
%% Author(s): Lorenzo A. Ricciardi (2015-2016)
% email: lorenzo.ricciardi@strath.ac.uk
%
%% References


% Run MACS
[memory,nfeval,ener]=macs7v16OC(fitnessfcn,[],LB,UB,options,[],[],varargin{1:end});

% Output
x    = memory(:,1:length(LB));
fval = memory(:,length(LB)+1:end-2);

output.memory               = memory;
output.nfeval               = nfeval;
output.ener                 = ener;
exitflag                    = 1;

end
