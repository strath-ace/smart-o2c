function [memories_out, memories, archivebest, population_evolution, vval_evolution, B_mean, delta_local, inite, iglob, options, exitflag] = ...
         MP_AIDEA(fname, vlb, vub, pop, options, varargin)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
% =========================================================================
% MP-AIDEA-ALR: Multi-Population Adaptive Inflationary Differential
% Evolution Algorithm with Adaptive Local Restart.
% - Adaptation of CR and F inside the single populations
% - Adaptation of the dimension of the bubble for the local restart
% - Adaptive local/global restart
% Reference: M. Di Carlo, M. Vasile, E. Minisci, "Multi-Population
% Adaptive Inflationary Differential Evolution Algorithm with Adaptive
% Local Restart", IEEE Congress on Evolutionary Computation, CEC 2015,
% Sendai, Japan
% =========================================================================
% =========================================================================
%    INPUT
%           fname   = function handle to cost function
%           vlb     = lower boundaries
%           vub     = upper boundaries
%           pop     = initial population matrix
%                     pop is a 3D matrix where the third dimension
%                     identifies the number of populations
%           options = structure containing options for the problem
%
%   OUTPUT
%          memories_out = archive containing the best individual, for each
%                         population, at each value of function evaluations
%                         chosen to record the results obtained
%          memories   = archive containing all the local minima plus the
%                       individuals of the 
%                       population at each restart. The solutions are sorted
%                       from the best to the worst. One cell for each
%                       population
%          archivebest = archive of local minima (obtained by local search)
%          population_evolution = initial and final population obtained at
%                                 each DE step for each population
%          vval_evolution = for each population, DE parameters 
%          B_mean   = mean value of the matrix for the adaptation of delta
%                     local
%          delta_local = value of delta_local for each population
%          inite    = number of local restart for each population
%          iglob    = number of global restart for each population
%          options
%          exitflag = exitflag from fmincon
% =========================================================================
% AIDEA:
% (c) Edmondo Minisci and Massimiliano Vasile 2013
% MP-AIDEA:
% (c)  Marilena Di Carlo 2015
% email: marilena.di-carlo@strath.ac.uk
% =========================================================================

%% Input check on algorithm parameters
% Check if the user defined all the required options for MP-AIDEA. Options
% not defined are assigner default parameters

% Check if number of population is compatible with defined input
if size(pop,3) == 1 && (~isfield(options,'delta_global') || isempty(options.delta_local))
    error('The adaptation of delta_local is not possible using only one population. Increase the number of populations or define a value for options.delta_local')
end

if size(pop,3) == 1 && (~isfield(options,'max_LR') || isempty(options.max_LR))
    error('The adaptation of the local restart is not possible using only one population. Increase the number of populations or define a value for options.max_LR')
end

% Check if delta_global was provided by the user
if ~isfield(options,'delta_global') || isempty(options.delta_global) 
    warning('delta_global (options.delta_global) is not defined or is empty. Default value of 0.1 is used. Press any key to continue or stop MP-AIDEA (Ctrl+c) and define a value for options.delta_global');
    pause
    options.delta_global = 0.1;
end

% Check if convergence threshold was provided by the user
if ~isfield(options,'rho') || isempty(options.rho)
    warning('rho (options.rho) is not defined by the user or is empty. The default value 0.2 is used. Press any key to continue or stop MP-AIDEA (Ctrl+c) and define a value for options.rho');
    pause
    options.rho = 0.2;
end

% Chek if user choose DE strategy
if ~isfield(options,'DE_strategy')  || isempty(options.DE_strategy)  
    warning('The DE strategy (options.DE_strategy) is not defined by the user. The default value is used (DE/Rand and DE/CurrentToBest). Press any key to continue or stop MP-AIDEA (Ctr+c) and define options.DE_strategy');
    pause
    options.DE_strategy = 1;

end

% Chek if user choose probability of DE strategy
if ~isfield(options,'prob_DE_strategy') || isempty(options.prob_DE_strategy)
    warning('The probablity of DE strategies (options.prob_DE_strategy) is not defined by the user. The default value 0.5 is used. Press any key to continue or stop MP-AIDEA (Ctrl+c) and define options.prob_DE_strategy');
    pause
    options.prob_DE_strategy = 0.5;
end

% Chek if user choose CRF
if ~isfield(options,'dd_CRF') || isempty(options.dd_CRF)
    warning('options.dd_CRF is not defined by the user. The default value 3 is used. Press any key to continue');
    pause
    options.dd_CRF = 3;

end

% Check on maximum number of function evaluations
if ~isfield(options,'nFeValMax') && isempty(options.nFeValMax)
    warning('The maximum number of function evaluations (options.nFeValMax) is not defined by the user. A default value of 100000 is used. Press any key to continue or stop MP-AIDEA (Ctrl+c) and define options.nFeValMax');
    pause
    options.nFeValMax = 100000;
end


% check on text flag
if ~isfield(options,'text') || isempty(options.text)
    warning('The flag for text display (options.text) is not defined by the user. Text will not be displayed during the optimisation. Press any key to continue or stop MP-AIDEA and define a value for options.text');
    pause
    options.text = 0;
end

%% MP-AIDEA parameters and options

% Number of individuals of the populations
NP         = size(pop, 1);

% Neighborhood limit for GLOBAL restart
expstepglo = options.delta_global;

% Probability of using strategy best
P6         = options.prob_DE_strategy;

% Limit on delta f for CR update (3)
dd_limit   = options.dd_CRF;

% Maximum number of function evaluations - NEW
nFeValMax  = options.nFeValMax;

% DE strategy
DE_strategy = options.DE_strategy;

% Problem dimension
D          = length(vlb);

% Number of populations
pop_number = size(pop,3);

% Size of the search space
DELTA      = vub-vlb;

% Used for the global restart; is the minimal distance from the centers of the clusters
distrmin   = sqrt(D*(expstepglo)^2);

% Initial cost function minimum value (to be compared to cost function
% values corresponding to population elements)
fmin       = 1e15;


%% Initialization - clear here

% Number of local restart performed by each population
inite      = zeros(1,pop_number);

% Number of global restart performed by each population
iglob      = zeros(1,pop_number);

% Initialize matrix containing all the local minima found by all the
% populations
archiveALL = [];

% Initialize matrix containing the mean value of the matrix B (adaptation
% of delta_local)
B_mean = [];

% Initialize matrices for the collection of the best members of the
% populaiton and the corresponding local minima and for the collection of
% the dimension of the region of attraction of each local minimum (used to
% perform adaptive local search and local/global restart)
ArchiveBM_LM = [];
d_region_attraction = [];

% Val is a variable with a number of row equal to the number of elements in
% the population and a number of column equal to the number of populations.
% For each population and each element of the population, it collect the
% cost function value
Val     = zeros(NP, pop_number);

% Best value of the function
BestVal = zeros(1, pop_number);

% Best member of the population
BestMem  = zeros(pop_number, D);

% Number of function evaluation for each population
nFeVal  = zeros(1, pop_number);

% Variable for the definition of contraction of each population. Contracted
% population are identified by 1 while non-contracted population are 0.
contraction = zeros(1, pop_number);

% Number of DE step performed by each population - re-initialized when
% contraction conditions are reached by all the population
pop_step = zeros(1, pop_number);

% Number of not-NaN elements in the matrix archivebest for each population.
% Values are summed also after population re-initialization.
archivebest_elements = zeros(1,pop_number);

% Vector for the identification of the populations that are waiting for the
% others population to be contracted. Each column is 1 when the
% corresponding population is waiting.
waiting = zeros(1,pop_number);


% Index for the creation of the matrix for the adaptation of the bubble
% dimension for the local restart.
% Initialize matrix_bubble_local to zero so that after the first local
% seach of each population, the matrix for the adaptation of delta_local
% can be defined
matrix_bubble_local  = 0;
% Initialize matrix_bubble_global to zero so that after the first global
% restart of each population, the matrix for the adaptation of delta_local
% can be defined
matrix_bubble_global = 0;

% ??
var1 = zeros(1,pop_number);
var2 = zeros(1,pop_number);
var3 = zeros(1,pop_number);


% Initialise vector
for i_pop_number = 1 : pop_number
    
    % Variable with number of cells equal to the number of the populations.
    % Each cells saves the individual of the population (but not their
    % function values) for each population, before and after each DE_step
    % (a DE step is until contraction of the population). The saved individual
    % will be therefore those at the begin of the DE and those at the end
    % of the process, at contraction
    population_evolution{i_pop_number} = [];
    
    % For each population, vval_evolution stores information during the DE
    % process. Information are not saved only after contraction of the
    % population, but all along the DE evolution (for each parents and
    % offspring generation). Each row of vval_evolution saves:
    % 1st column: DE iteration number
    % 2nd column: minimum value of objective function for the current
    % individuals of the current generation of the considered population
    % 3rd column: mean value of the objective function 
    % 4th column: maximum value of the objective function
    % 5th column: parameter contraction
    % 6th column: parameter contraction
    vval_evolution{i_pop_number}       = [];
    
    % For each population, saves the individuals and the objective function
    % value of each individual. Differ from population_evolution for the
    % presence of the objective function value and because it includes also
    % the local minima (not only the individuals of the DE population).
    % Moreover its values are sorted from the minimum to the maximum
    % objective function value (losing therefore the sequentiality of the
    % individuals)
    memories{i_pop_number}             = [];
    
    % Archive of local minima found by local search 
    archivebest{i_pop_number}          = [];
    
    % Value of delta_local for each population
    delta_local{i_pop_number}          = [];
end

% Solutions will be saved when fraction of the maximum number of function evaluations will be
% reached. This fractions are defined in input.record. An example could be:
% record = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
record = options.record;

% Initialise step for the saving of solution at fraction of maximum number
% of function evaluations
step = 1;


%% Main loop
% Exit when the maximum number of function evaluations has been reached:
while sum(nFeVal) < nFeValMax
    
    %% DE loop for all the populations 
    % Exit from this cycle when population contracts (break condition) or when
    % the maximum number of function evaluations is reached
    while sum(nFeVal) < nFeValMax
        
        % Cycle over all the populations
        for i_pop_number = 1 : pop_number
            
            % =================================================================
            % If one of the population has already performed a global restart,
            % it will wait for all the other population to perform it before
            % being allowed to enter the following cycle.
            % This condition is expressed through the following "if" condition.
            % iglob is a variable with a number of columns equal to the number
            % of populations which gives, for each population, the number of
            % time that that population has performed a global restart.
            if iglob(1,i_pop_number) == min(iglob)
                
                
                % =============================================================
                % PARENTS --> OFFSPRING
                % =============================================================

                % Advance from parents to offspring
                if options.text 
                    disp('---------------------------------------------------------------------------------------')
                    disp('DIFFERENTIAL EVOLUTION')
                    disp(['Population number: ' num2str(i_pop_number)]);
                    disp('---------------------------------------------------------------------------------------')
                end

                % Save current population 
                population_evolution{i_pop_number} = [population_evolution{i_pop_number}; pop(:,:,i_pop_number)]; 
                
                % ---------------------------------------------------------
                % Run DE until contraction
                % ---------------------------------------------------------
                [BestMem(i_pop_number,:),... % Best member of the current population
                 BestVal(1,i_pop_number),... % Function value associated to the best member of the current population
                 nFeVal,...                  % Function evaluations
                 pop(:,:,i_pop_number),...   % Current population after DE
                 Val(:,i_pop_number),...     % Function value for the individual of the current population
                 ~,...                       % number of generations for contraction (not used)
                 vval_DE,...                 % Vector with useful quantities for the evolution of the population
                 new_elements ...            % Number of individuals considered for the current population before DE stopped because nFeValMax was reached
                 exitflag ...                % Flag for fmincon - used if cycle exit before fmincon               
                 ] = DE_step(fname, ...                 % Handle to function to optimise
                             vlb,...                    % Lower boundaries
                             vub,...                    % Upper boundaries
                             pop(:,:,i_pop_number),...  % Current population
                             nFeVal,...                 % Vector with number of function evaluations for each population
                             i_pop_number,...           % index of current population
                             options,...                 % options (containing: rho, P6, dd, DE_strategy)
                             nFeValMax,...              % Maximum number of functions evaluations
                             varargin{:});
                

                      
                % Save population and other parameters of the population
                population_evolution{i_pop_number} = [population_evolution{i_pop_number}; pop(:,:,i_pop_number)];  
                vval_evolution{i_pop_number}       = [vval_evolution{i_pop_number};       vval_DE];
                
                % ---------------------------------------------------------
                % Check condition related to maximum number of function
                % evaluations. 
                % Is it necessary to exit the cycle?
                % ---------------------------------------------------------
                % If number of function evaluation is greater than number
                % of function evaluation for which it is necessary to save
                % results:
                if sum(nFeVal) >= record(step) * nFeValMax
                    
                    % If the total maximum number of function evaluations
                    % has been reached
                    if sum(nFeVal)>= nFeValMax
                        
                        % For cycle over all the populations
                        for i_population = 1 : pop_number
                            
                            % Update populations with index lower than "i_pop_number" 
                            % (current population) with all the elements (NP)
                            if i_population < i_pop_number
                                
                                % Add all the individuals of the population
                                memories_update = [pop(:,:,i_population) Val(:,i_population)];
                                
                                % Add memories_update
                                memories{i_population} = [memories{i_population}; memories_update];
                                
                                % Sort individuals based on their objective
                                % function value
                                memories{i_population} = sortrows(memories{i_population}, D+1);
                                
                                % For the current step of the record
                                % process, and for the current population,
                                % take the best individual of the current
                                % population (1st row, individuals were
                                % sorted)
                                memories_out(step, :, i_population) = memories{i_population}(1,:);
                                
                                % Update the population i_pop_number with a number of
                                % elements equal to new_elements (and no more than that)
                            elseif i_population == i_pop_number
                                
                                % What is "new_elements"? If the function
                                % DE_step was interrupted because the
                                % maximum number of function ev was
                                % reached, not all the individuals of the
                                % population could have been evaluated.
                                % "new_elements" is the number of evaluated
                                % individual before DE_step was interrupted
                                memories_update = [pop(1 : new_elements, : ,i_population) Val(1:new_elements,i_population)];
                                
                                % Update memories
                                memories{i_population} = [memories{i_population}; memories_update];
                                
                                % Sort row based on objective function
                                % value
                                memories{i_population} = sortrows(memories{i_population}, D+1);
                                
                                memories_out(step, :, i_population) = memories{i_population}(1,:);
                                
                                % The remaining population did not have the opportunity to do a DE
                            elseif i_population > i_pop_number
                                
                                if ~isempty( memories{i_population} )
                                    
                                    memories{i_population} = sortrows(memories{i_population}, D+1);
                                    memories_out(step, :, i_population) = memories{i_population}(1,:);
                                    
                                end
                                
                            end
                            
                            % end of for cycle over number of population
                        end
                        
                        % Since this is the case in which the total maximum
                        % number of function evaluation is reached, exit
                        % the function:
                        return
                        
                        % else, if total maximum number of function 
                        % evaluations has not been reached, but it is only
                        % necessary to save results for record:
                    else

                        for i_population = 1 : i_pop_number
                            
                            memories_update = [pop(:,:,i_population) Val(:,i_population)];
                            memories{i_population} = [memories{i_population}; memories_update];
                            memories{i_population} = sortrows(memories{i_population}, D+1);
                            memories_out(step, :, i_population) = memories{i_population}(1,:);
                        end
                        
                        for i_population = i_pop_number + 1 : pop_number
                            
                            
                            if ~isempty(memories{i_population})
                                memories{i_population} = sortrows(memories{i_population}, D+1);
                                memories_out(step, :, i_population) = memories{i_population}(1,:);
                            end
                        end
                        
                        % Increase step of record by one
                        step = step + 1;
                    end
                end
                % ---------------------------------------------------------------------
                % End of check condition related to maximum number of function evaluations
                % ---------------------------------------------------------------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                
             
                % The current population has finished the DE cycle above
                % and has therefore reached contraction. Put to 1 the flag
                % to show that this has happened
                contraction(1,i_pop_number ) = 1;
                % =================================================================
                % end of -----if iglob(1,i_pop_number) == min(iglob)
            end
            % =================================================================
            
            % =====================================================================
            % end of ---------for i_pop_number = 1 : pop_number
        end
        % ====================================================================
        
        
        % =====================================================================
        % Identification of the condition for the exit of this cycle
        % We exit when all the population have reached contraction condition
        % Pay attention that when one or more population have already been
        % global restarted (and are waiting for the other to be globally
        % restared too) we do not have to have 1 in each element of the vector
        % contraction because we will have 0 in correspondence of the
        % population already globally restarted 
        % =====================================================================     
        
        % 1st case - all the population have had the same number of global
        % restart
        if rem(sum(iglob), pop_number) == 0
            
            if sum(contraction) == length(contraction)
                % all population are contracted
                
                % Bring the contraction values back to zero for the next cycle
                % and exit this cycle (go to local-global restart phase)
                contraction = zeros(1,pop_number);
                pop_step    = zeros(1,pop_number);
                break
            end
            
            % 2nd case - each population has had a different number of global
            % restart
        else
            
            % Compute how many population have been globally restarted up to
            % now
            number_global_pop = rem(sum(iglob), pop_number);
            
            % Exit the cycle if the number of contracted population is equal to
            % the number of population we are still considering in this cycle
            % that is, (pop_number - number_global_pop)
            if sum(contraction) >= ( pop_number - number_global_pop )
                
                % All the not-globally reinitialized population are contracted
                % Bring the contraction values back to zero for the next cycle
                % and exit this cycle (go to local-global restart phase)
                contraction = zeros(1,pop_number);
                pop_step    = zeros(1,pop_number);
                break
            end
        end
        
    end
    
    
    
    %% Local search and local restart / Global restart
    
    % =========================================================================
    % Runs local search and re-initialize population (locally or globally)
    % =========================================================================
    
    % Number of function evaluations available for the local search
    % algorithm    
    nFeValLS = round((nFeValMax - sum(nFeVal)) /pop_number);
   
    % If there are more than 1 population, use mechanism to adapt
    % delta_local and local/global restart
    if pop_number > 1 && isempty(options.max_LR)
        
        for i_pop_number = 1 : pop_number
            
            % Do not try to search for local minimum if the population has
            % been already globally restarted and is waiting for the other
            % populations to be globally restarted too
            if iglob(1,i_pop_number) == min(iglob)
                
                
                % Create a 3D matrix that collects all best members BM and local minima LM
                % ArchiveBM_LM = [BM f(BM); LM f(LM)]
                % for each best member from which a local seach is started and
                % its corresponding local minimum. Collect all of them in 3d
                % arrray
                if isempty(ArchiveBM_LM)
                    ArchiveBM_LM(1,:,1) = [BestMem(i_pop_number,:) BestVal(1,i_pop_number)];
                    ArchiveBM_LM(2,:,1) = zeros(1, D+1);
                else
                    ArchiveBM_LM = cat(3, ArchiveBM_LM, NaN*ones(size(ArchiveBM_LM,1), size(ArchiveBM_LM,2)) );
                    ArchiveBM_LM(1,:,end) = [BestMem(i_pop_number,:) BestVal(1,i_pop_number)];
                end
                
                % Dimension of the region of attraction of a local minimum
                d_region_attraction = cat(2,d_region_attraction,[NaN; NaN]);
                
                % Particular cases for local/global restart. See below
                case1 = ~exist('first_local_restart','var');
                case2 = exist('first_local_restart','var') && var1(1,i_pop_number) == 0;
                case3 = rem(sum(iglob),pop_number) == 0 && var2(1,i_pop_number) == 0 && any(iglob);
                case4 = rem(sum(iglob),pop_number) == 0 && var3(1,i_pop_number) == 0 && any(iglob);
                
                % CASE 1
                % If this is the first possibility to run a local search i.e.
                % no local search was ever conducted before: perform local
                % search and assume that the found local minimum is not inside
                % any region of attraction of already detected local minimum
                % (there are no stored local minimum). Do this for all the
                % populations
                % CASE 2
                % The local restart has been performed once for all the
                % populations and this is the second local search for the
                % considered population: realise local search, assume that we
                % are not in the basin of attraction and put var1 to 1
                % for the current population so that we can not enter this
                % case2 never again. The reason for this is that
                % CASE 3
                % All the populations have had the same number of global
                % restart (different from zero)
                % CASE 4
                % All the populations have had the same number of global
                % restart (different from zero)
                if case1
                    local_search = 1;
                    inside_region_attraction = 0;
                elseif case2
                    local_search = 1;
                    inside_region_attraction = 0;
                    var1(1,i_pop_number) = 1;
                elseif case3
                    local_search = 1;
                    inside_region_attraction = 0;
                    var2(1,i_pop_number) = 1;
                elseif case4
                    local_search = 1;
                    inside_region_attraction = 0;
                    var3(1,i_pop_number) = 1;
                else
                    
                    % If we are in none of the cases above, we can proceed as
                    % usual, by checking if the best member is inside the
                    % region of attraction of already detected local minima
                    
                    
                    % Very high distance from best member to local minimum (we
                    % are going to compute a minimum value below)
                    distance_BM_wrt_LM = 10^(34);
                    
                    % Compute minimum distance between current best member and previous
                    % local minima.
                    % The current best member is ArchiveBM_LM(1,1:D,end) while
                    % the previous local minimum is ArchiveBM_LM(2,1:3,i)
                    for i = 1 : size(ArchiveBM_LM,3) - 1
                        
                        % Minimum distance
                        distance_BM_wrt_LM = min(distance_BM_wrt_LM, norm(ArchiveBM_LM(1,1:D,end) - ArchiveBM_LM(2,1:D,i)));
                        
                        % As soon as the best member is within the region of attraction
                        % of some local minimum (d_region of attraction computed with
                        % 4 best members! not only one)
                        % The 3rd condition of the following line means that the
                        % current f(BM) has to be higher than the considered
                        % f(LM)
                        if distance_BM_wrt_LM <= d_region_attraction(1,i)  && ...
                                d_region_attraction(2,i) >= 4   && ...
                                ArchiveBM_LM(1,D+1,end) >= ArchiveBM_LM(2,D+1,i)
                            
                            % We are inside the region of attraction of a local
                            % minimum! Break the for cycle because there is no need to
                            % Check all the other minima
                            inside_region_attraction = 1;
                            
                            % Index to identify the minimum in whom basin of
                            % attraction we are
                            index_region = i;
                            break
                        else
                            inside_region_attraction = 0;
                        end
                        
                    end
                    
                end
                
                
                % If best member BM is inside the region of attraction of an
                % already detected local minimum:
                if inside_region_attraction
                    
                    % No local search:
                    warning_LS(i_pop_number) = 1;
                    
                    
                    % if the best member BM is NOT inside the region of attraction of an already
                    % detected local minimum or if for any of the cases it is
                    % necessary to perform a local search:
                elseif ~inside_region_attraction || local_search == 1
                    
                    
                    % ==================================================================
                    % Local Search
                    % ==================================================================
                    
                    % Number of function evaluations for the local search for
                    % the current population
                    nfev(1,i_pop_number) =round( max([min([300*D  nFeValLS]) D]));
                    
                    % fmincon finds a constrained minimum of a function of several
                    % variables.
                    % The start point is the minimum function value defined through BestMem.
                    % The solution is searched in the interval lb to ub.
                    % The otuputs are:
                    % - xgrad      location of the minimum
                    % - fvalgrad   value of the function at the minimum point
                    % - exitflag   describes the exit condition of fmincon (integer
                    %              identifying the reason the algorithm terminated)
                    % - output     structure with information about the optimization
                    %              (output.funcCount = number of function evaluations)
                     
                    foptionsNLP = optimset('Display','off','MaxFunEvals',nfev(1,i_pop_number),'LargeScale','off','FinDiffType','central','Algorithm','sqp');
                    
                    % Local search with fmincon
                    if ~isfield(options,'no_bounds')
                        [xgrad,fvalgrad,exitflag,output] = fmincon(fname,BestMem(i_pop_number,:),[],[],[],[],vlb,vub,[],foptionsNLP,varargin{:});
                    elseif isfield(options,'no_bounds') && options.no_bounds
                        % For problem with no bounds use fminunc instead
                        % than fmincon
                        [xgrad,fvalgrad,exitflag,output] = fminunc(fname,BestMem(i_pop_number,:),foptionsNLP,varargin{:});
                    end
                    
                    % Update number of function evaluations
                    nFeVal(1,i_pop_number) = nFeVal(1,i_pop_number) + output.funcCount;
                    
                    % If value of minimum find by the local search is lower
                    % than current best value:
                    if fvalgrad < BestVal(1,i_pop_number)
                        % substitute best value with function value of minimum
                        % found by fmincon
                        BestVal(1,i_pop_number) = fvalgrad;
                        % Substitute best member components with minimum found
                        % by fmincon
                        BestMem(i_pop_number,:) = xgrad;
                    end
                    if options.text
                        disp('---------------------------------------------------------------------------------------')
                        disp('LOCAL SEARCH')
                        disp(['Population number: ' num2str(i_pop_number)]);
                        disp('---------------------------------------------------------------------------------------')
                        disp(['Minimum: ' num2str(BestVal(1,i_pop_number))]);
                        disp(['Number of function evaluations (so far): ' num2str(sum(nFeVal))]);
                        disp('------------------------------------------------')
                    end
                    
                    % ==================================================================
                    % Archiving of solutions
                    % ==================================================================
                    % Add local minimum to the archive
                    ArchiveBM_LM(2,:,end) = [xgrad fvalgrad];
                    
                    % isnan is true for Not-A-Number
                    if isnan(BestVal(1,i_pop_number))
                        error('NaN value in the objective function')
                    else
                        % archivebest is updated at each local minimum search with the new
                        % bestvalue and the number of function evaluation (only if the
                        % minimum thus found gives a value of the function lower than
                        % the one given by the best element of the population)
                        
                        archivebest_update = [BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                        
                        archivebest{i_pop_number} = [archivebest{i_pop_number}; archivebest_update];
                        
                        % Archivebest collects values separately for each
                        % population. ArchiveALL collects all the minimum
                        % together (used by clustering mean shift for
                        % clustering of local minima and global restart of the
                        % population)
                        archiveALL = [archiveALL; BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                        
                        % We have added a new element to the archivebest matrix for the
                        % i_pop_number population
                        archivebest_elements(1,i_pop_number) = archivebest_elements(1,i_pop_number) + 1;
                        
                        % The population used until this point is going to be
                        % re-initialized. Before this happens, save the population into
                        % the matrix memories - therefore memories contains the
                        % population at the moment of the contraction of the population
                        memories_update = [pop(:,:,i_pop_number) Val(:,i_pop_number);...
                            xgrad                 fvalgrad];
                        
                        memories{i_pop_number} = [memories{i_pop_number}; memories_update];
                        
                    end
                    
                    
                    % ==================================================================
                    % Update Best
                    % ==================================================================
                    
                    % If the new BestVal value found for the current contracted
                    % population is not better than the previous ones (fmin) the
                    % population keeps being re-initialized around the previous xmin
                    % value
                    
                    if BestVal(1,i_pop_number) < fmin
                        
                        fmin = BestVal(1,i_pop_number);
                        xmin = BestMem(i_pop_number,:);
                        
                    end
                    
                    % What happens to the following line if BestVal is not below the fmin
                    % value and therefore xmin is not defined? This does not happen - see
                    % the value of fmin!!!
                    % Reference point for the re-initialization of the
                    % population
                    xref(i_pop_number,:) = xmin;
                    
                    % Dimension of the region of attraction: distance between
                    % best member (ArchiveBM_LM(1,1:D,end)) and local minimum
                    % (ArchiveBM_LM(2,1:D,end)) for the current best member and
                    % local minimum (identified by "end")
                    d_region_attraction(:,end) = [norm(ArchiveBM_LM(1,1:D,end)-ArchiveBM_LM(2,1:D,end)); 1];
                    
                    %
                    warning_LS(i_pop_number) = 0;
                    
                    % Check if the detected minimum was already detected by
                    % looking in the archive of minima
                    for i = 1 : size(ArchiveBM_LM,3)-1
                        
                        % Compute distance of current minimum from all the
                        % minima in the archive of minima detected by the local
                        % search
                        distance_LM_wrt_LM   = norm( ArchiveBM_LM(2,1:D,end) - ArchiveBM_LM(2,1:D,i) );
                        
                        % Compute the distance in objective function between
                        % the current minimum and all the minima in the archive
                        difference_LM_wrt_LM = norm( ArchiveBM_LM(2,D+1,end) - ArchiveBM_LM(2,D+1,i) );
                        
                        % If the distance between local minimum is lower than a
                        % certain value and if the difference in objective
                        % function is lower than the difference in objective
                        % function between current  best member and best member
                        % of the considered local minimum, the local minimum
                        % was already detected!
                        if (  distance_LM_wrt_LM   <= norm(DELTA) * 0.001 ) && ...
                                ( difference_LM_wrt_LM <= norm(ArchiveBM_LM(1,D+1,end)-ArchiveBM_LM(1,D+1,i)))
                            
                            % Update dimension of the region of attraction of
                            % current minimum
                            d_region_attraction(1,end) = min(d_region_attraction(1,end), d_region_attraction(1,i));
                            
                            % Update dimension of the region of attraction of
                            % minimum "i"
                            d_region_attraction(1,i) = d_region_attraction(1,end);
                            
                            % Update number of times that minimum was detected
                            % both for minimum in "i" and current minimum "end"
                            d_region_attraction(2,i)   = d_region_attraction(2,i)  + 1;
                            d_region_attraction(2,end) = d_region_attraction(2,end) + 1;
                            
                        end
                        
                        % end for cycle over best members and local minimum in archive
                    end
                    
                    % end of condition "inside region of attraction"
                end
                
                % end of if to check if population is waiting for global
                % restart
            end
            
            
            % ---------------------------------------------------------------------
            % Check condition related to maximum number of function evaluations
            % ---------------------------------------------------------------------
            if sum(nFeVal) >= record(step) * nFeValMax
                
                for i_pop_number = 1 : pop_number
                    
                    if ~isempty(memories{i_pop_number})
                        memories{i_pop_number}    = sortrows(memories{i_pop_number},D+1);
                        memories_out(step,:,i_pop_number) = memories{i_pop_number}(1,:);
                    end
                end
                
                if sum(nFeVal)>= nFeValMax
                    return
                else
                    step = step + 1;
                end
            end
            
            % End of cycle for number of populations
        end
        

        
        
        
    else
       
        % IF the population is only one or if max_LR is an input, run local search always before going into
        % local/gobal restart
        
        for i_pop_number = 1 : pop_number
            
            % Number of function evaluations for the local search for
            % the current population
            nfev(1,i_pop_number) =round( max([min([300*D  nFeValLS]) D]));
            
            foptionsNLP = optimset('Display','off','MaxFunEvals',nfev(1,i_pop_number),'LargeScale','off','FinDiffType','central','Algorithm','sqp');
            
            % Local search with fmincon
            [xgrad,fvalgrad,exitflag,output] = fmincon(fname,BestMem(i_pop_number,:),[],[],[],[],vlb,vub,[],foptionsNLP,varargin{:});
            
            % Update number of function evaluations
            nFeVal(1,i_pop_number) = nFeVal(1,i_pop_number) + output.funcCount;
            
            % If value of minimum find by the local search is lower
            % than current best value:
            if fvalgrad < BestVal(1,i_pop_number)
                % substitute best value with function value of minimum
                % found by fmincon
                BestVal(1,i_pop_number) = fvalgrad;
                % Substitute best member components with minimum found
                % by fmincon
                BestMem(i_pop_number,:) = xgrad;
            end
            if options.text
                disp('---------------------------------------------------------------------------------------')
                disp('LOCAL SEARCH')
                disp(['Population number: ' num2str(i_pop_number)]);
                disp('---------------------------------------------------------------------------------------')
                disp(['Minimum: ' num2str(BestVal(1,i_pop_number))]);
                disp(['Number of function evaluations (so far): ' num2str(sum(nFeVal))]);
                disp('------------------------------------------------')
            end
            
            % ==================================================================
            % Archiving of solutions
            % ==================================================================
            % Add local minimum to the archive
            ArchiveBM_LM(2,:,end) = [xgrad fvalgrad];
            
            % isnan is true for Not-A-Number
            if isnan(BestVal(1,i_pop_number))
                error('NaN value in the objective function')
            else
                % archivebest is updated at each local minimum search with the new
                % bestvalue and the number of function evaluation (only if the
                % minimum thus found gives a value of the function lower than
                % the one given by the best element of the population)
                
                archivebest_update = [BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                
                archivebest{i_pop_number} = [archivebest{i_pop_number}; archivebest_update];
                
                % Archivebest collects values separately for each
                % population. ArchiveALL collects all the minimum
                % together (used by clustering mean shift for
                % clustering of local minima and global restart of the
                % population)
                archiveALL = [archiveALL; BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                
                % The population used until this point is going to be
                % re-initialized. Before this happens, save the population into
                % the matrix memories - therefore memories contains the
                % population at the moment of the contraction of the population
                memories_update = [pop(:,:,i_pop_number) Val(:,i_pop_number);...
                    xgrad                 fvalgrad];
                
                memories{i_pop_number} = [memories{i_pop_number}; memories_update];
                
            end
            
            
            % ==================================================================
            % Update Best
            % ==================================================================
            
            % If the new BestVal value found for the current contracted
            % population is not better than the previous ones (fmin) the
            % population keeps being re-initialized around the previous xmin
            % value
            
            if BestVal(1,i_pop_number) < fmin
                
                fmin = BestVal(1,i_pop_number);
                xmin = BestMem(i_pop_number,:);
                
            end
            
            % What happens to the following line if BestVal is not below the fmin
            % value and therefore xmin is not defined? This does not happen - see
            % the value of fmin!!!
            % Reference point for the re-initialization of the
            % population
            xref(i_pop_number,:) = xmin;
            
            
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If the dimension for the local restart is not defined by the user,
        % create matrix to adapt it and adapt it when required:
        if isempty(options.delta_local)
            
            % =========================================================================
            % Creation of a matrix for the adaptation of the dimension of the bubble
            % =========================================================================
            % This matrix has to be created considering the local minima found by the
            % population in two occasion:
            % - at the first local search, before the first local restart
            % - after each global restart of all the populations
            % The code for now separates these two events but the aim is to merge the
            % two while and for loop that follows in a single one.
            
            while matrix_bubble_local == 0
                
                
                % Distance between local minima
                min_minima_distance = 10^60;
                max_minima_distance = 0;
                
                % After the first local search of each population, the number of
                % local minima will be equal to the number of populations
                for i = 1 : pop_number-1
                    
                    for k = i+1 : pop_number
                        
                        % Distance between local minima
                        minima_distance = norm(BestMem(i,:) - BestMem(k,:));
                        
                        % Minimum distance between local minima
                        min_minima_distance = min(min_minima_distance, minima_distance);
                        
                        % Maximum distance between local minima
                        max_minima_distance = max(max_minima_distance, minima_distance);
                        
                    end
                    
                end
                
                %
                max_minima_distance = max_minima_distance + norm(DELTA)*0.1;
                
                % INITIALIZATION OF THE MATRIX FOR THE ADAPTATION OF THE DIMENSION OF
                % THE BUBBLE
                % The dimension of the bubble for the local restart is defined in the
                % interval min_minima_distance to max_minima_distance
                
                % If D = 3 (dimension of the problem) B will be (D+1)x3 where the second
                % and third coulmns correspond to the parameter p for the adaptation of the
                % dimension of the bubble and to the iteration number
                delta     = (max_minima_distance - min_minima_distance)/(D/2);
                B_1st     = (min_minima_distance : delta : max_minima_distance)';
                B         = [B_1st zeros(size(B_1st)) zeros(size(B_1st))];
                
                % The matrix for the dimension of the bubble for the first local
                % restart has been created. This cycle will not be entered anymore.
                matrix_bubble_local = 1;
                
                % Shall we adapt the matrix right now after this cycle? No, if we are
                % going to have the first local restart.
                first_local_restart = 1;
                
            end
            
            % This matrix has to be created also (using different condition) after each
            % global restart of ALL POPULATIONS, considering all the elements in archivebest.
            
            % If all the population have had the same number of global restart
            % (that is different from zero)
            if rem(sum(iglob), pop_number) == 0  && any(iglob)  && matrix_bubble_global == 0
                
                % Distance between all the local minima
                for i = 1 : size(archiveALL,1)
                    
                    for k = i+1 : size(archiveALL,1)
                        
                        % Distance between local minima
                        minima_distance = norm(archiveALL(i,:) - archiveALL(k,:));
                        
                        % Minimum and maximum distance between local minima
                        min_minima_distance = min(min_minima_distance, minima_distance);
                        max_minima_distance = max(max_minima_distance, minima_distance);
                        
                    end
                    
                end
                max_minima_distance = max_minima_distance + norm(DELTA)*0.1;
                
                
                % INITIALIZATION OF THE MATRIX FOR THE ADAPTATION OF THE DIMENSION OF
                % THE BUBBLE
                %         % The dimension of the bubble for the local restart is defined in the
                % interval min_minima_distance to max_minima_distance
                
                % If D = 3 (dimension of the problem) B will be (D+1)x3 where the second
                % and third coulmns correspond to the parameter p for the adaptation of the
                % dimension of the bubble and to the iteration number
                delta     = (max_minima_distance - min_minima_distance)/(D/2);
                B_1st     = (min_minima_distance : delta : max_minima_distance)';
                B         = [B_1st zeros(size(B_1st)) zeros(size(B_1st))];
                
                % We need to perform this matrix generation only once. Therefore bring
                % the following value to 1 after the matrix has been generated.
                matrix_bubble_global = 1;
                
                % Shall we adapt the matrix right now after this cycle? No, if we are
                % going to have the first local restart.
                first_local_restart = 1;
                
                
                
            end
            
            
            % =========================================================================
            % Adaptation of the dimension of the bubble
            % =========================================================================
            % The adaptation is performed before each local restart but not AT the
            % first local restart, when we have still no values for the comparison. At
            % the same time, after the global restart, the matrix is not adapted but
            % we wait for the subsequent local restart to adapt it
            
            % No matrix adaptation if we are going to have the first local restart
            if first_local_restart == 0
                
                for i_pop_number = 1 : pop_number
                    
                    % The adaptation is performed only for the population which are not
                    % already globally restarted and are waiting for other populations
                    % to globally restart too
                    if iglob(1,i_pop_number) == min(iglob)
                        
                        % Difference between previous local minimum and current local
                        % minimum - this has to be different from zero for a bubble to
                        % be sufficiently big
                        p = norm(archivebest{i_pop_number}(end,1:D) - archivebest{i_pop_number}(end-1,1:D) );
                        
                        for ic = 1 : size(B,1)
                            
                            % If the previous difference between minima is smaller than
                            % the new one
                            if B(ic,2) < p
                                
                                % Substitue the previous value of the bubble with the one used
                                % for the considered population
                                B(ic,1) = bubble_dimension(1,i_pop_number);
                                
                                % Associate the corresponding p and iter value
                                B(ic,2) = p;
                                B(ic,3) = pop_step(1,i_pop_number);
                                
                                break
                            end
                            
                        end
                        
                    end
                    
                    B = sortrows(B, [2 3]);
                    
                    % Save value of B mean for output
                    B_mean = [B_mean; mean(B(:,1))/norm(DELTA)];
                end
                
            end
            
            % BUBBLE DIMENSIONS SAMPLED FROM PARZEN DISTRIBUTION
            % Bv is a pop_number*1 matrix which, for each population elements (row)
            % contains the bubble dimension value
            Bv = parzenself_k([], B(:,1), ones(size(B(:,1))), pop_number, 'norm', 0);
            % Avoid Bv outside some boundaries?!?!
            
            % End of condition if isempty(options.delta_local)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % =========================================================================
    % Local and Global Restart
    % =========================================================================
    
    for i_pop_number = 1 : pop_number
        
        % If value of delta_local is not provided by the user:
        if isempty(options.delta_local)
            % Link each bubble dimension value to a population
            bubble_dimension(1, i_pop_number) = Bv(i_pop_number,1);
        else
            % If the dimension of the bubble is an input of the user:
            bubble_dimension(1, i_pop_number) = options.delta_local * norm(DELTA);
        end
        
        % save dimension of the bubble for the local restart for each
        % population
        delta_local{i_pop_number} = [delta_local{i_pop_number}; bubble_dimension(1,i_pop_number)/norm(DELTA)];

        
        % Do not try to perform local or global restart if the population has
        % been already globally restarted and is waiting for the other
        % populations to be globally restarted too
        if iglob(1,i_pop_number) == min(iglob)
            
            % Perform local restart when:
            % - after 1st local search of all the populations (warning_LS
            % is not defined)
            % - at any point, when warning_LS ~= 1, meaning that the best
            % member of the population is not inside the basin of
            % attraction of any local minimum, and therefore local search
            % and local restart are realized
            % - the population is only one and its number of local restart
            % is below the maximum number

            if (pop_number>1 && isempty(options.max_LR) && ~exist('warning_LS','var') ) || ...
               (exist('warning_LS','var') && warning_LS(i_pop_number) ~= 1) || ...
               (pop_number == 1 && inite(1,i_pop_number) < options.max_LR) || ...
               (pop_number>1 && ~isempty(options.max_LR) && inite(1,i_pop_number) < options.max_LR)
                
                
                % --------- LOCAL RESTART
                
                % EACH population is restarted in a local bubble around xref.
                % The boundaries of the bubble are defined throught the following
                % two lines
                
                % For problem with no bounds
                if isfield(options,'no_bounds') && options.no_bounds
                    
                    % Create an hypercube defined by the minimum and maximum value
                    % of the elements during the evolution, so far
                    max_NoBounds = -1e16;
                    min_NoBounds =  1e16;
                    
                    for index_temp = 1 : pop_number
                        max_NoBounds = max(max_NoBounds, max(max(memories{index_temp}(:,1:D))));
                        min_NoBounds = min(min_NoBounds, min(min(memories{index_temp}(:,1:D))));
                    end
                    
                    vub_NoBounds = max_NoBounds * ones(1,D);
                    vlb_NoBounds = min_NoBounds * ones(1,D);
                    
                    XVminl = max([xref(i_pop_number,:) - bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vlb_NoBounds]);
                    XVmaxl = min([xref(i_pop_number,:) + bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vub_NoBounds]);
                    
                else
                    % For problem with bounds
                    
                    XVminl = max([xref(i_pop_number,:) - bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vlb]);
                    XVmaxl = min([xref(i_pop_number,:) + bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vub]);
                end
         
                % Population re-initialized in the bubble
                pop(:,:,i_pop_number) = lhsdesign(NP,D,'criterion','maximin').*repmat(XVmaxl-XVminl,NP,1)+repmat(XVminl,NP,1);
                
                % Increase number of local restart for current population
                % by 1
                inite(1,i_pop_number) = inite(1,i_pop_number) + 1;
                
                if options.text
                    disp('---------------------------------------------------------------------------------------')
                    disp('LOCAL RESTART')
                    disp(['Population number: ' num2str(i_pop_number)]);
                    disp(['Number of local restart for current population: ' num2str(inite(1,i_pop_number))]);
                    disp('---------------------------------------------------------------------------------------')
                end
                
                % if warning_LS=1, no local search and local restart is
                % performed, but rather population is globally restarted
                % Global restart is performed also when population is only
                % 1 and the maximum number of local restart has been
                % reached for that population
            elseif (exist('warning_LS','var') && warning_LS(i_pop_number) == 1) || ...
                   (pop_number == 1  && inite(1,i_pop_number) >= options.max_LR) || ... 
                   (pop_number > 1  && ~isempty(options.max_LR) && inite(1,i_pop_number) >= options.max_LR) 
                
                % --------- GLOBAL RESTART
                
                % Count the number of global restart for each population
                iglob(1,i_pop_number) = iglob(1,i_pop_number) + 1;
                
                % 
                if rem(sum(iglob),pop_number) == 0
                    var2 = zeros(1,pop_number);
                    var3 = zeros(1,pop_number);
                end
                
                % The current population is going to be globally restarted.
                % From now on, it will wait for all the other populations to be
                % globally restarted before starting another DE evolution.
                % Therefore the corresponding value goes to 1 in the vector
                % waiting.
                waiting(1,i_pop_number) = 1;
                
                % Clustering of the archived local minima
                [ClusterCenters,DataClusters,DatainClusters] = Clustering_MeanShift(archiveALL(:,1:D)',distrmin);
               
                % Number of clusters
                usc = size(ClusterCenters,2);
                
                % Center of each cluster
                meanclu=ClusterCenters';
                
                pop0=[];
                
                while size(pop0,1) < NP
                                        
                    % For problem with no bounds
                    if isfield(options,'no_bounds') && options.no_bounds
                        xref0 = lhsdesign(NP,D,'criterion','maximin').* repmat(vub_NoBounds-vlb_NoBounds,NP,1) + repmat(vlb_NoBounds,NP,1);                
                    else
                        xref0 = lhsdesign(NP,D,'criterion','maximin').* repmat(vub-vlb,NP,1) + repmat(vlb,NP,1);
                    end
    
                    for i = 1 : NP
                        
                        dmin=1e15;
                        
                        % Checking distance condition for each member of the
                        % population xref0 from each cluster (for cycle going from
                        % 1 to the number of clusters usc). For the considered i-th
                        % member of the population, the minumum distance is saved
                        % in dmin
                        for j = 1 : usc
                            dmin = min([dmin norm(meanclu(j,:)-xref0(i,:))]);
                        end
                        
                        
                        % If the minimum distance of the considered element of the
                        % population with respect to the cluster centers is below a
                        % certian value the element is accepted, otherwise is
                        % eliminated.
                        % distrmin is the minimal distance from the centers of the clusters
                        
                        % A condition has been added to the following IF cycle
                        % to take into account the possibility to overcome the
                        % number of elements of the population - rivedere e
                        % dire a max
                        if dmin >= distrmin  && (size(pop0,1)<NP)
                            pop0=[pop0; xref0(i,:)];
                        else
                            %                         disp('population element deleted')
                        end
                        
                    end
                    
                end
                
                pop(:,:,i_pop_number) = pop0;
                
                
                % Current population contraction status is back to zero
                contraction(1,i_pop_number) = 0;
                
                % Before starting with the next series of local restart bring
                % the fmin values back to their original values
                fmin=1e15 * ones(1,pop_number);
                
                if pop_number == 1
                   % If one population is considered, put inite to 0 after
                   % global restart, in order to be able to compare the
                   % number of local restart to options.max_LR
                    inite(1,i_pop_number) = 0;
                    
                end
                
                
                if options.text
                    disp('---------------------------------------------------------------------------------------')
                    disp('GLOBAL RESTART')
                    disp(['Population number: ' num2str(i_pop_number)]);
                    disp(['Number of global restart for current population: ' num2str(iglob(1,i_pop_number))]);
                    disp('---------------------------------------------------------------------------------------')
                end
            end
            
        end
        % ---------------------------------------------------------------------
        

        % We need to exit the part of the code that performs local and global
        % restart when all population are globally restarted.
        % We cannot use the parameter iglob for this check.
        % We use waiting, which is the vector which identifies with 1 the
        % population that are currently waiting for all the others to be globally
        % restarted.
        % When all the population are in a waiting conditions, we can exit this
        % part of the code, but waiting has to be brought back to zero before.
        if sum(waiting) == pop_number
            waiting = zeros(1, pop_number);
            matrix_bubble_global = 0;
            break
        end
        
        % ========================================================================
        % end of -------- for i_pop_number = 1 ; pop_number
    end
    % =========================================================================
        
    
    % We have had the first local restart, so at the next cycle we can adapt
    % the matrix for the dimension of the bubble
    first_local_restart = 0;
   
    
    % ========================================================================
    % end of -------- while 1
    % that comprise both the DE step + contraction condition evaluation and the
    % local/global restart part
end
% =========================================================================


% =========================================================================
% end of the function
end
% =========================================================================






%% Differential Evolution with adaptive F and CR 

function [BestMem, BestVal, nFeVal, pop, Val, iter, vval, new_elements, exitflag] = DE_step(fname, ...
         lb, ub, pop, nFeVal, i_pop_number, options, nFeValMax, varargin)

%    INPUT
%           fname        : function handle to cost function
%           lb,ub        : lower and upper boundaries
%           pop          : input population (parents)
%           nFeVal       : vector with the number of functions evaluations for
%                          each population
%           i_pop_number : number of the current population
%           options      : structure containing:
%                          mmdist       : convergence threshold (for population contraction)
%                          P6           : probability of running selected DE strategies
%                          dd_limit     : limit value of dd for adaptation of CR and F
%                          DE_strategy  : integer to identify selected DE strategies    
%                          CR           : CR for non-adaptive run
%                          F            : F for non-adaptive run
%           nFeValMax    : maximum number of function evaluations allowed
%           varargin     : additional inputs
%
%   OUTPUT
%           BestMem       : best solution vector
%           BestVal       : best solution value
%           nFeVal        : number of function evaluations
%           pop           : population (children at contraction of the
%                           population)
%           Val           : cost function associated to pop
%           iter          : number of generations for contraction
%           vval          : vector with useful information about the
%                           evolution of the population
%           new_elements  :



% (c) Edmondo Minisci and Massimiliano Vasile 2013
%     Marilena Di Carlo, 2015

% =========================================================================
% Initialization
% =========================================================================
% Contraction threshold
mmdist = options.rho;

% Probability of using strategy best
P6         = options.prob_DE_strategy;

% Limit on delta f for CR update (3)
dd_limit   = options.dd_CRF;

% DE strategy
DE_strategy = options.DE_strategy;

% CR and F for the DE defined by the user. If empty, they will be adapted
CR_user = options.CR;
F_user  = options.F;

% Number of new individuals evaluated before reaching the maximum number of
% function evaluation
new_elements = 0;

% 
exitflag = 0;

% The population is composed by NP elements with dimension D
[NP,D]    = size(pop);

% Val is what in idea2.m was called fitness - vector collecting the
% objective function value for all the individuals of the population
Val       = zeros(1,NP);


% Mean element of the population
MeanPop   = mean(pop);

% Maximum distance between elements of the population and their mean value
mmdistm = -1.e10;

for imm = 1 : size(pop)
    dista = norm(pop(imm,:) - MeanPop);
    if dista > mmdistm
        mmdistm = dista;
    end
end

% Maximum distance before the offspring generation
mmdistm0 = mmdistm;

% Row line summarizing the function values
% mmdistm/mmdistm0 is equal to 1 for the parent generation and becomes
% lower than 1 (and always smaller) during the contraction of the
% population
vval=[0 min(Val) mean(Val) max(Val) mmdistm mmdistm/mmdistm0];


% =========================================================================
% Function Evaluation for Initial Population (Parents)
% =========================================================================
% start with first population member
ibest   = 1;
Val(1)  = feval(fname, pop(ibest,:), varargin{:});
BestVal = Val(1);                 % best objective function value so far
nFeVal(1,i_pop_number)  = nFeVal(1,i_pop_number) + 1;

if sum(nFeVal) >= nFeValMax
    
    iter = 0;
    new_elements = 1;
    
    BestMem = pop(ibest,:);
    BestVal = Val(1);

    return
end

for i = 2 : NP                        % check the remaining members
    Val(i) = feval(fname, pop(i,:), varargin{:});
    nFeVal(1,i_pop_number) = nFeVal(1,i_pop_number) + 1;
    
    if (Val(i) < BestVal)           % if member is better
        ibest   = i;                 % save its location
        BestVal = Val(i);
    end
    
    if sum(nFeVal) >= nFeValMax
        new_elements = i;
%         vval = 0;
        iter = 0;
        
        BestMem = pop(ibest,:);

        return
    end
    
end

CurrBest = pop(ibest,:);            % best member of current iteration
BestMem = CurrBest;                 % best member ever


% =========================================================================
% DE Initialization
% =========================================================================

RotIndArr = (0:1:NP-1);               % rotating index array (size NP)

iter = 1;

nostop = 1;


%% Initialization of CRF to uniform distribution

% Initialize CRF only if CR and F are to be adapted (that is, if they are
% not defined)
if isempty(CR_user) && isempty(F_user)
    
    % Crossover probability CR is defined in the interval 0.1 to 0.99
    % CRa is a matrix with ((D/2)+1) rows and 3 columns
    % If D = 3 CRa will be
    % CRa = [0.1  0  0;...
    %        CR2  0  0;...
    %        CR3  0  0;...
    %        0.99 0  0];
    delta = (.99-.1)/(D/2);
    CRa_1st   = (0.1:delta:0.99)';          % 1st column of matrix CRa
    % CRa_1st = linspace(0.1, 0.99, fix(sqrt(D + 1)))';
    CRa       = [CRa_1st zeros(size(CRa_1st)) zeros(size(CRa_1st))];
    
    % Differential weight F is defined in the interval -1 to 1
    % Fa is a (D/2+1)*3 matrix
    delta     = (1-(-1))/(D/2);
    Fa_1st    = (-1:delta:1)';             % 1st column of matrix Fa
    % Fa_1st = linspace(-1, 1, fix(sqrt(D + 1)))';
    Fa        = [Fa_1st zeros(size(Fa_1st)) zeros(size(Fa_1st))];
    
    CRFa = [];
    
    % Construnction of a regular mesh in the space with CR between 0.1 and 0.99
    % and F between -1 and 1. The third and foruth coulmns correspond to the
    % parameter dd and to the iteration number
    % At the end of the following for cycle the matrix CRFa will be (if D=3)
    % CRFa = [CR1 F1 0 0;...
    %         CR1 F2 0 0;...
    %         CR1 F3 0 0;...
    %         CR1 F4 0 0;...
    %         CR2 F1 0 0;...
    %         CR2 F2 0 0:...
    %         ......
    %         CR3 F1 0 0;...
    %         .....
    %         CR4 F4 0 0];
    % and therefore its dimensions will be [(D/2+1)x(D/2+1)]*4
    for im = 1 : size(CRa(:,1))
        for in = 1 : size(Fa(:,1))
            CRFa=[CRFa; CRa(im,1) Fa(in,1) 0 0];
        end
    end
    
    
end

%% Main DE loop
while nostop
    
    % If CR and F are not defined by the user
    if isempty(CR_user) && isempty(F_user)
        
        % =====================================================================
        % CR and F are sampled from the Parzen distribution
        % =====================================================================
        % CRFv is a NP*2 matrix which, for each population elements (row)
        % contains the CR (1st column) and F (2nd column) values
        CRFv = parzenself_k([], CRFa(:,1:2), ones(size(CRFa(:,1))), NP, 'norm', 0);
        
        % Separate CR from F values
        % At the end of the following lines:
        % CRv -> NP*1
        % Fv  -> NP*1
        CRv = CRFv(:,1);
        Fv  = CRFv(:,2);
        
        % ---------------------------------------------------------------------
        % Avoid CR outside the boundaries 0.1 to 0.99
        % ---------------------------------------------------------------------
        % If some of the obtained CRv associated to a population element is
        % lower than 0.1, its values is brought back to 0.1
        % At the same time, is some element have a value bigger than 0.99, it
        % is brought back to 0.99
        % The matrix CRv as defined in the following lines is indeed
        % CRv = [CR1 CR2 CR3 ..... CRn;...
        %        0.1 0.1 0.1       0.1]
        % where CR1, CR2 etc represent the CR values associated to the 1st, 2nd
        % and so on individual of the population.
        % The max function returns a row vector containing the maximum value in
        % each column
        % The same holds for the min comparison
        % At the end of the following two lines CRv will be a NP column vector
        % (it is transposed back at the end of the lines!)
        CRv = (max ( [CRv'; 0.1  * ones(size(CRv')) ] ) )';
        CRv = (min ( [CRv'; 0.99 * ones(size(CRv')) ] ) )';
        % CR is of size [D,NP]
        CR  = repmat(CRv, 1, D);
        
        Fv  = (max ( [Fv'; -1*ones(size(Fv')) ] ) )';
        Fv  = (min ( [Fv';  1*ones(size(Fv')) ] ) )';
        % F is of size [D, NP]
        F   = repmat(Fv, 1, D);
        
        
        % If CR and F are defined by the user
    else
        
        
        CR = repmat(CR_user, NP, D);
        F  = repmat(F_user, NP, D);
    end
    
    
    
    % =====================================================================
    % Generation of the offspring
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Probabilistic selection of the strategy 
    % ---------------------------------------------------------------------
    if rand(1,1) < P6
        strategy = 1;
    else
        strategy = 2;
    end
    
    
    % ---------------------------------------------------------------------    
    % Creation of the offspring population
    % --------------------------------------------------------------------- 
    
    % Save old population in "popold"
    popold = pop;                 
    
    % Index pointer array
    IndPointArr = randperm(4);              
    
    % Shuffle locations of vectors
    IndArr1  = randperm(NP);             
    
    % Rotate indices by IndPointArr(1) positions
    RotIndArr2 = rem(RotIndArr + IndPointArr(1), NP); 
    
    IndArr2  = IndArr1(RotIndArr2+1);               
    RotIndArr2 = rem(RotIndArr+IndPointArr(2),NP);
    IndArr3  = IndArr2(RotIndArr2+1);
   
    PopMat1 = popold(IndArr1,:);             % shuffled population 1
    PopMat2 = popold(IndArr2,:);             % shuffled population 2
    PopMat3 = popold(IndArr3,:);             % shuffled population 3
    
    % Population filled with the best member of the last iteration
    % ---Old implementation:
%     for i = 1 : NP                      
%         BestMat(i,:) = CurrBest;         
%     end
    % ---New implementation
    BestMat = repmat(CurrBest, NP, 1);
    

    % Mask e
    % To each individual of the population a different value of
    % CR is associated!
    MaskInterPop = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise
    MaskOldPop = MaskInterPop < 0.5;         % inverse mask to MaskInterPop
    
   
    % ---------------------------------------------------------------------
    % DE_strategy = 1 : DE/CurrentToBest, DE/Rand
    % ---------------------------------------------------------------------
    
    if DE_strategy == 1
        
        if (strategy == 1) % DE/CurrentToBest
            InterPop = popold  + F.*(BestMat - popold) + F.*(PopMat1 - PopMat2);                      % differential variation
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;                                   % crossover
        elseif (strategy == 2) % DE/Rand
            InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
        end
        
        
        % ---------------------------------------------------------------------
        % DE_strategy = 2 : DE/Best, DE/Rand
        % ---------------------------------------------------------------------
    elseif DE_strategy == 2
        
        if (strategy == 1) % DE/Best
            InterPop = BestMat + F.*(PopMat1 - PopMat2);                        % differential variation
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;             % crossover
        elseif (strategy == 2)   % DE/Rand                                                     
            InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
        end
        
        
    end
    
    % =====================================================================
    % Select which vectors are allowed to enter the new population
    %======================================================================
    if isempty(CR_user) && isempty(F_user)
        CRFa = sortrows(CRFa, [3 4]);
    end
    
    % When I exit this function I want to return the maximum dd value for
    % all the elements of the population. Therefore the part of the code
    % related to the compuation of dd has been changed wrt to the original
    % aidea code. In particular, dd is initialized to zero and at each
    % step of the for cycle that analyze the element of the population, the
    % maximum value of dd is returned.
    % The aim is to return the maximum dd value to perform the adaptation
    % of the CRF matrix.
    
    
    for i = 1 : NP
        
        % =================================================================
        % Every individual violating the boundaries is projected back into
        % the boundaries
        % =================================================================
        
        % Distance of the i-th element of the population from the lower
        % boundaries
        dInterPopl = InterPop(i,:) - lb;
        
        % How many component of the i-th element of the population are
        % below the lower boundary?
        ElTooLow  = find(dInterPopl < 0);
        lElTooLow = length(ElTooLow);
        
        % Distance of the i-th element of the population from the upper
        % boundaries
        dInterPopr = ub - InterPop(i,:);
        
        % How many component of the i-th element of the population are
        % above the lower boundary?
        ElTooHigh  = find(dInterPopr < 0);
        lElTooHigh = length(ElTooHigh);
        
        % Bring back the considered element into the boundaries of the
        % problem:
        
        % If problem is not without bounds, bring individuals that are
        % outside of the boundaries, back into the search space
        if ~isfield(options,'no_bounds')
            
            % -----------------------------------------------------------------
            % Random re-initialisation
            % -----------------------------------------------------------------
            %         if lElTooLow > 0
            %            InterPop(i,ElTooLow) = lb(ElTooLow) + rand(1,lElTooLow).*(ub(ElTooLow) - lb(ElTooLow));
            %         end
            %
            %         if lElTooHigh > 0
            %            InterPop(i,ElTooHigh) = lb(ElTooHigh) + rand(1,lElTooHigh).*(ub(ElTooHigh) - lb(ElTooHigh));
            %         end
            
            % -----------------------------------------------------------------
            % Bounce back re-initialisation
            % -----------------------------------------------------------------
            if lElTooLow > 0
                InterPop(i,ElTooLow) = ( popold(i,ElTooLow) + lb(ElTooLow) ) / 2;
            end
            
            if lElTooHigh > 0
                InterPop(i,ElTooHigh) = ( popold(i,ElTooHigh) + ub(ElTooHigh) ) / 2;
            end
            
        elseif isfield(options,'no_bounds') && options.no_bounds
            % If problem is without bounds
        end
        
        % =================================================================
        % Value of the function f for each new element of the population
        % =================================================================
        
        TempVal = feval(fname, InterPop(i,:), varargin{:});   
        
        % Increase number of function evalutations
        nFeVal(1,i_pop_number)  = nFeVal(1,i_pop_number) + 1;
        
        
        if (TempVal <= Val(i))  % if competitor is better than parent
            
            dd = abs((TempVal-Val(i)));
            
            % Replace parent vector with child (for new iteration)
            pop(i,:) = InterPop(i,:);
            
            % Save value
            Val(i) = TempVal;
            
            if isempty(CR_user) && isempty(F_user)
                % -------------------------------------------------------------
                % Adaptation of CR and F
                % -------------------------------------------------------------
                for ic = 1 : size(CRFa,1)
                    
                    % If the previous difference between parent and offspring
                    % function value is smaller than the new one (under the
                    % condition that the offspring perform better than the
                    % parent - we are already in an if condition, that is if
                    % TempVal <= Val(i) )
                    if CRFa(ic,3) < dd
                        
                        % The substitution of CR is subject to another
                        % limitation
                        if abs(dd) > dd_limit
                            
                            % Substitution of CR
                            CRFa(ic,1) = CRv(i);
                        end
                        
                        % Substitue the previous value of F with the one used
                        % for the considered offspring
                        CRFa(ic,2) = Fv(i);
                        
                        % Associate the corresponding dd and iter value
                        CRFa(ic,3) = dd;
                        CRFa(ic,4) = iter;
                        break
                    end
                end
                
            end
            
            % If the value of the function for the considered child
            % individual is better than best value ever..
            if (TempVal < BestVal)     
                BestVal = TempVal;            % new best value
                BestMem = InterPop(i,:);      % new best parameter vector ever
            end
        end
        
        if sum(nFeVal) >= nFeValMax
            new_elements = i;
            
            return
        end
        
    end %---end for imember=1:NP
    
    
    CurrBest = BestMem;       % freeze the best member of this iteration for the coming
    % iteration. This is needed for some of the strategies.
    
    
    % =====================================================================
    % Evaluate maximum distance between elements of the population
    % =====================================================================
    
    % Mean element of the population
    MeanPop = mean(pop);
    
    mmdistm = -1.e10;
    
    for imm = 1 : size(pop)
        
        % Distance of individual from mean individual
        dista = norm(pop(imm,:) - MeanPop);
        if dista > mmdistm
            mmdistm = dista;
        end
    end
    
    % mmdistm is the maximum distance between population mean and elements
    % of the population

    % Collects results of the considered offspring
    vvalr = [iter min(Val) mean(Val) max(Val) mmdistm mmdistm/mmdistm0];
    vval = [vval; vvalr];
    
    % Increase number of generations by one
    iter = iter + 1;
    
    % Exit when one of the following conditions happens:
    % - the population contracts (mmdistm/max(vval(:,5))<mmdist)
    % - the maximum number of function evaluations has been reached
    % (sum(nFeVal)>nFeValMax)
    % - the DE is running from more than (D/10)*100 generations    
    nostop = mmdistm/max(vval(:,5))>mmdist  && sum(nFeVal)<nFeValMax && iter < (D/10)*100 ;

    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.text
    disp(['Contraction at DE generation: ' num2str(iter)]);
    disp(['Minimum: ' num2str(BestVal)]);
    disp('------------------------------------------------')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



