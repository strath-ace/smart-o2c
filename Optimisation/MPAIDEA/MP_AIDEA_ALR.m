function [memories_out,vval, B_mean,bubble, archivebest,options, exitflag] = MP_AIDEA_ALR(fname,vlb, vub, pop, options, varargin)
%
% =========================================================================
% - Multipopulation Adaptive AIDEA
% - DE strategy randomly chosen between DE/best and DE/rand, or DE/best and
%   DE/arch, or DE/arch and DE/rand
% - Adaptation of CR and F inside the single populations
% - Adaptation of the dimension of the bubble
% - DE strategy randomly chosen between DE/best and DE/rand
% =========================================================================
% =========================================================================
%
%   [memories,archivebest,options] = MP_AIDEA_SINGLE_AdaptBubble(fname,pop,vlb,vub,options,varargin)
%
%    INPUT
%           fname   = function handle to cost function
%           pop     = initial population matrix
%                     if pop is an empty or partial matrix then aidea
%                     will generate the missing elements
%                     pop is a 3D matrix where the third dimension identify
%                     the numbers of populations
%           vlb     = lower boundaries
%           vub     = upper boundaries
%           options = options vector
%                     options(4)  = number of elements of the population
%                     options(12) = neighborhood limit for GLOBAL restart
%                     options(27) = population contraction threshold
%                     options(31) = limit number of local restart
%                     options(34) = probability of using strategy
%                                   DE/best/1/bin (rather than
%                                   DE/rand/1/bin) or others
%                     options(35) = limit on delta f for CR update (CRC)
%                     options(36) = 1 to get a video, 0 otherwise
%
%   OUTPUT
%          memories   = archive containing all the local minima plus the
%                       population for each restart. The solutions are sorted
%                       from the best to the worst
%          archivebest = archive of local minima

% (c) Edmondo Minisci and Massimiliano Vasile 2013
%     Marilena Di Carlo 2014
%
% =========================================================================
B_mean = [];
bubble.pop1 = [];
bubble.pop2 = [];
bubble.pop3 = [];
bubble.pop4 = [];

%% Default parameters

% Check if delta_global was provided by the user
if ~isfield(options,'delta_global')    
    warning('\delta_{global} (options.delta_global) not defined - Default value is used');
    options.delta_global = max(vub-vlb)/50;
end

% Check if convergence threshold was provided by the user
if ~isfield(options,'rho')
    warning('\rho (options.rho) not defined by the user - Default value is used');
    options.rho = max(vub-vlb)/100;
end

% Chek if user choose DE strategy
if ~isfield(options,'DE_strategy')    
    warning('DE strategy (options.prob_DE_strategy) not defined by the user - Default value is used');
    options.prob_DE_strategy = 0.5;

end

% Chek if user choose probability of DE strategy
if ~isfield(options,'prob_DE_strategy')
    warning('Probablity of DE strategy (options.prob_DE_strategy) not defined by the user - Default value is used');
    options.prob_DE_strategy = 0.5;
end

% Chek if user choose CRF
if ~isfield(options,'dd_CRF')
    warning('options.dd_CRF not defined by the user - Default value is used');
    options.dd_CRF = 3;

end

% Check on maximum number of function evaluations
if ~isfield(options,'nFeValMax')
    warning('Maximum number of function evaluations (options.nFeValMax) not defined by the user - Default value is used');
    options.nFeValMax = 100000;
end

% check on plot flag
if ~isfield(options,'plot_flag')
        warning('Flag for plots (options.plot_flag) not defined by the user - Default value is used');
    options.plot_flag = 0;

end

% check on plot flag
if ~isfield(options,'text')
        warning('Flag for text (options.text) not defined by the user - Default value is used');
    options.text = 0;
end

%% MP-AIDEA parameters and options

% Number of individuals of the populations
NP         = size(pop, 1);

% Neighborhood limit for GLOBAL restart
expstepglo = options.delta_global;

% Convergence threshold
mmdist     = options.rho;

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


%% Initialization

% Number of local restart performed by each population
inite      = zeros(1,pop_number);

% Number of global restart performed by each population
iglob      = zeros(1,pop_number);

% Initialize archivebest (matrix for the collection of all the local minima
% found using fmincon) and memories (matrix for the collection of all
% population elements and local minima)
archivebest=[];
memories   =[];

archiveALL = [];

% Initialize matrices for the collection of the best members of the
% populaiton and the corresponding local minima and for the collection of
% the dimension of the region of attraction of each local minimum
ArchiveBM_LM = [];
d_region_attraction = [];

% Val is a variable with a number of row equal to the number of elements in
% the population and a number of column equal to the number of populations.
% For each population and each element of the population, it collect the
% cost function value
Val     = zeros(NP, pop_number);

% Best value of the function
BestVal = zeros(1,pop_number);

% Number of function evaluation for each population
nFeVal  = zeros(1,pop_number);

% Best member of the population
BestMem  = zeros(pop_number,D);

% Maximum distance between population elements (for contraction
% definition)
mmdistm = -1.e10*ones(1,pop_number);

%
vval = zeros(2,6,pop_number);
vval_new = zeros(2,6,pop_number);

% Variable for the definition of contraction of each population. Contracted
% population are characterized by 1 while non-contracted population are 0.
contraction = zeros(1,pop_number);

% Number of DE step performed by each population - re-initialized when
% contraction conditions are reached by all the population
pop_step = zeros(1,pop_number);

% Number of not-NaN elements in the matrix archivebest for each population.
% Values are summed also after population re-initialization.
archivebest_elements = zeros(1,pop_number);

% Parameter for the identification of the number of the rows in the matrix
% archivebest related to each new series of local optimizer (see code
% below)
row_number=0;

% Vector for the identification of the populations that are waiting for the
% others population to be contracted. Each column is 1 when the
% corresponding population is waiting.
waiting = zeros(1,pop_number);


% Index for the creation of the matrix for the adaptation of the bubble
% dimension for the local restart
matrix_bubble_local  = 0;
matrix_bubble_global = 0;
minima = [];

var1 = zeros(1,pop_number);
var2 = zeros(1,pop_number);
var3 = zeros(1,pop_number);



for i_pop_number = 1 : pop_number
    vval_NEW{i_pop_number} = [];
end


step = 1;
record = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];



%% Main loop
% Exit when the maximum number of function evaluations has been reached
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
            % being allowed to enter the following cycle and therefore advance
            % in the DE algorith.
            % This condition is expressed through the following line.
            % iglob is a variable with a number of columns equal to the number
            % of populations which gives, for each population, the number of
            % time that that population has performed a global restart.
            if iglob(1,i_pop_number) == min(iglob)
                
                
                % =============================================================
                % PARENTS --> OFFSPRING
                % =============================================================
                
                % Advance from parents to offspring
                if options.text == 1
                    disp('---------------------------------------------------------------------------------------')
                    disp('DIFFERENTIAL EVOLUTION')
                    disp(['Population number: ' num2str(i_pop_number)]);
                    disp('---------------------------------------------------------------------------------------')
                end

                
                [BestMem(i_pop_number,:),... % Best member of the current population
                 BestVal(1,i_pop_number),... % Function value associated to the best member of the current population
                 nFeVal,...                  % Function evaluations
                 pop(:,:,i_pop_number),...   % Current population after DE
                 Val(:,i_pop_number),...     % Function value for the individual of the current population
                 iter,...                    % 
                 vval_DE,...                 % 
                 new_elements ...            %
                 ] = DE_step(fname, ...                 % Handle to function to optimise
                             vlb,...                    % Lower boundaries
                             vub,...                    % Upper boundaries
                             pop(:,:,i_pop_number),...  % Current population
                             nFeVal,...                 % Vector with number of function evaluations for each population
                             i_pop_number,...           % index of current population
                             mmdist,...                 % Contraction threshold
                             P6,...                     % Probability of using the two selected DE strategies
                             archivebest,...            % 
                             dd_limit,...               % CRF
                             DE_strategy,...            % Selected DE strategy 
                             nFeValMax,...              % Maximum number of functions evaluations
                             varargin{:});
                
                vval_NEW{i_pop_number} = [vval_NEW{i_pop_number}; vval_DE];
                
                % ---------------------------------------------------------------------
                % Check condition related to maximum number of function evaluations
                % ---------------------------------------------------------------------
                if sum(nFeVal) >= record(step)*nFeValMax
                    
                    if sum(nFeVal)>= nFeValMax
                        for i_population = 1 : pop_number
                            if i_population < i_pop_number
                                % Update each memories with NP new elements
                                memories(end+1:end+NP,:,1:pop_number)= NaN*ones(NP,D+1,pop_number);
                                memories_update = [pop(:,:,i_population) Val(:,i_population)];
                                memories(end-NP+1:end,:,i_population)=memories_update;
                                memories(:,:,i_population)    = sortrows(memories(:,:,i_population),D+1);
                                memories_out(step,:,i_population) = memories(1,:,i_population);
                            elseif i_population == i_pop_number
                                % Update this populatio with new_elements
                                % new elements
                                memories(end+1:end+new_elements,:,1:pop_number)= NaN*ones(new_elements,D+1,pop_number);
                                memories_update = [pop(1:new_elements,:,i_population) Val(1:new_elements,i_population)];
                                memories(end-new_elements+1:end,:,i_population)=memories_update;
                                memories(:,:,i_population)    = sortrows(memories(:,:,i_population),D+1);
                                memories_out(step,:,i_population) = memories(1,:,i_population);
                            elseif i_population > i_pop_number
                                % Take best element
                                memories(:,:,i_population)    = sortrows(memories(:,:,i_population),D+1);
                                memories_out(step,:,i_population) = memories(1,:,i_population);
                            end
                        end
                        sum(nFeVal)
                        return
                    else
                        memories_out(step,:,:) = NaN * ones(1,D+1,pop_number);
                        for i_population = 1 : i_pop_number
                            memories(end+1:end+NP,:,1:pop_number)= NaN*ones(NP,D+1,pop_number);
                            memories_update = [pop(:,:,i_population) Val(:,i_population)];
                            memories(end-NP+1:end,:,i_population)=memories_update;
                            memories(:,:,i_population)    = sortrows(memories(:,:,i_population),D+1);
                            memories_out(step,:,i_population) = memories(1,:,i_population);
                        end
                        
                        for i_population = i_pop_number + 1 : pop_number
                            memories(:,:,i_population)    = sortrows(memories(:,:,i_population),D+1);
                            memories_out(step,:,i_population) = memories(1,:,i_population);
                        end
                        step = step + 1;
                    end
                end
                
                
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
        % population already globally contracted - this could be made in a
        % better way, for example not updating contraction to zero when exiting
        % from this cycle and going to the local-global restart cycle
        % and bringing it back to zero only if the local restart is performed.
        % NOOOOOO perche il DE step si realizza solo se contraction e` uguale a
        % zero, quindi lasciarlo cosi com'e` e spiegare che non si puo fare
        % altrimenti perche se no contraction non e` zero
        % =====================================================================
        % the following part can be realized in a better way - maybe without
        % distinguishing between 1st and 2nd case... there should be a way to
        % merge the two cases into one
        
        
        % 1st case - all the population have had the same number of global
        % restart
        if rem(sum(iglob),pop_number) == 0
            
            if sum(contraction) == length(contraction)
                % all population are contracted
                % Bring the contraction values back to zero for the next cycle
                % and exit this cycle (go to local-global restart phase)
                contraction = zeros(1,pop_number);
                pop_step = zeros(1,pop_number);
                break
            end
            
            % 2nd case - each population has had a different number of global
            % restart
        else
            
            % Compute how many population have been globally restarted up to
            % now
            number_global_pop = rem(sum(iglob),pop_number);
            
            % Exit the cycle if the number of contracted population is equal to
            % the number of population we are still considering in this cycle
            % that is, (pop_number - number_global_pop)
            if sum(contraction) >= ( pop_number - number_global_pop )
                
                % All the not-globally reinitialized population are contracted
                % Bring the contraction values back to zero for the next cycle
                % and exit this cycle (go to local-global restart phase)
                contraction = zeros(1,pop_number);
                pop_step = zeros(1,pop_number);
                break
            end
        end
        
    end
    
    % Once all the population are contracted, vval is redefined in order to be
    % ready for a new session of DE step after the restart of all the
    % populations.
    % The values used of vval are to be related only to a single series of DE
    % steps, and the values are not to be kept to avoid problems at the next
    % series of DE step.
    vval = zeros(2,6,pop_number);
    
    
    %% Local search and local restart / Global restart
    
    % =========================================================================
    % Runs local search and re-initialize population
    % Aggiungere condizione che dopo un certo numero di local restart fa il
    % global
    % =========================================================================
    
    % A line to the archivebest matrix has to be added for each local restart
    % optimizer procedure in order to allow the storage of the new results.
    % However, this line has to be added only one time for each major step
    % (parallel local restart of all the population) in such a way as to have a
    % new line in each 3D matrix for each population.
    % The initial values are set to NaN so that if the considered population is
    % globally restarted and is waiting for the other population to globally
    % restart too, it will have a NaN value on the corresponding line and when
    % the sortrows is computed, all this values will go to the last rows of the
    % matrix
    
    for i_pop_number = 1 : pop_number
        
        if iglob(1,i_pop_number) == min(iglob)
            archivebest(end+1,1:D+2,1:pop_number) = NaN*ones(1,D+2,pop_number);
            memories(end+1:end+1+NP,:,1:pop_number)= NaN*ones(NP+1,D+1,pop_number);
            break
        end
        
    end
    
    nFeValLS = round((nFeValMax - sum(nFeVal)) /pop_number);
    
    for i_pop_number = 1 : pop_number
        
        
        
        % Do not try to search for local minimum if the population has
        % been already globally restarted and is waiting for the other
        % populations to be globally restarted too
        if iglob(1,i_pop_number) == min(iglob)
            
            
            % Create a matrix that collects all best members BM and local minima LM
            % ArchiveBM_LM = [BM f(BM); LM f(LM)]
            if isempty(ArchiveBM_LM)
                ArchiveBM_LM(1,:,1) = [BestMem(i_pop_number,:) BestVal(1,i_pop_number)];
                ArchiveBM_LM(2,:,1) = zeros(1, D+1);
            else
                ArchiveBM_LM = cat(3, ArchiveBM_LM, NaN*ones(size(ArchiveBM_LM,1), size(ArchiveBM_LM,2)) );
                ArchiveBM_LM(1,:,end) = [BestMem(i_pop_number,:) BestVal(1,i_pop_number)];
            end
            
            d_region_attraction = cat(2,d_region_attraction,[NaN; NaN]);
            
            
            case1 = ~exist('first_local_restart','var');
            case2 = exist('first_local_restart','var') && var1(1,i_pop_number) == 0;
            case3 = rem(sum(iglob),pop_number) == 0 && var2(1,i_pop_number) == 0 && any(iglob);
            case4 = rem(sum(iglob),pop_number) == 0 && var3(1,i_pop_number) == 0 && any(iglob);
            
            
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
                
                
                distance_BM_wrt_LM = 10^(34);
                
                % Compute minimum distance between current best member and previous
                % local minima
                for i = 1 : size(ArchiveBM_LM,3)-1
                    
                    distance_BM_wrt_LM = min(distance_BM_wrt_LM, norm(ArchiveBM_LM(1,1:D,end) - ArchiveBM_LM(2,1:D,i)));
                    
                    % As soon as the best member is within the region of attraction
                    % of some local minimum (d_region of attraction computed with
                    % two best members! not only one)
                    % The 3rd condition of the following line means that the
                    % current BM has to be higher than the considered LM
                    if distance_BM_wrt_LM <= d_region_attraction(1,i)  && d_region_attraction(2,i) >= 4   && ArchiveBM_LM(1,D+1,end) >= ArchiveBM_LM(2,D+1,i)
                        
                        inside_region_attraction = 1;
                        index_region = i;
                        break
                    else
                        inside_region_attraction = 0;
                    end
                    
                end
                
            end
            
            
            % Two things can happen if inside a region of attraction of a
            % previous detected local minimum: either I have already tried
            % escaping that region of attraction using the local restart or
            % I have not
            
            if inside_region_attraction
                
                % No local search
                warning_LS(i_pop_number) = 1;
                
                
                
            elseif ~inside_region_attraction || local_search == 1
                
                
                % ==================================================================
                % Local Optimizer
                % ==================================================================
                
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
                [xgrad,fvalgrad,exitflag,output] = fmincon(fname,BestMem(i_pop_number,:)',[],[],[],[],vlb,vub,[],foptionsNLP,varargin{:});
                
                nFeVal(1,i_pop_number) = nFeVal(1,i_pop_number) + output.funcCount;
                
                if fvalgrad < BestVal(1,i_pop_number)
                    BestVal(1,i_pop_number) = fvalgrad;
                    BestMem(i_pop_number,:) = xgrad;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             disp('---------------------------------------------------------------------------------------')
                %             disp('LOCAL SEARCH')
                %             disp(['Population number: ' num2str(i_pop_number)]);
                %             disp('---------------------------------------------------------------------------------------')
                %             pause(5)
                %             disp(['Minimum: ' num2str(BestVal(1,i_pop_number))]);
                %             disp('------------------------------------------------')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ==================================================================
                % Archiving of solutions
                % ==================================================================
                % Add local minimum to the archive
                ArchiveBM_LM(2,:,end) = [xgrad' fvalgrad];
                
                
                % isnan is true for Not-A-Number
                if isnan(BestVal(1,i_pop_number))
                    disp('NAN')
                    pop
                    BestMem
                    [isave nFeVal]
                    pause
                else
                    % archivebest is updated at each local minimum search with the new
                    % bestvalue and the number of function evaluation (only if the
                    % minimum thus found gives a value of the function lower than
                    % the one given by the best element of the population)
                    
                    archivebest_update = [BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                    archivebest(end,:,i_pop_number) = archivebest_update;
                    
                    archiveALL = [archiveALL; BestMem(i_pop_number,:) BestVal(1,i_pop_number) nFeVal(1,i_pop_number)];
                    
                    % We have added a new element to the archivebest matrix for the
                    % i_pop_number population
                    archivebest_elements(1,i_pop_number) = archivebest_elements(1,i_pop_number) + 1;
                    
                    % The population used until this point is going to be
                    % re-initialized. Before this happens, save the population into
                    % the matrix memories - therefore memories contains the
                    % population at the moment of the contraction of the population
                    memories_update = [pop(:,:,i_pop_number) Val(:,i_pop_number);...
                        xgrad' fvalgrad];
                    memories(end-NP:end,:,i_pop_number)=memories_update;
                end
                
                
                % ==================================================================
                % Update Best
                % ==================================================================
                
                % If the new BestVal value found for the current contracted
                % population is not better than the previous ones (fmin) the
                % population keeps being re-initialized aroung the previous xmin
                % value
                
                if BestVal(1,i_pop_number) < fmin
                    
                    dd = abs((BestVal(1,i_pop_number)-fmin)/fmin);
                    fmin = BestVal(1,i_pop_number);
                    xmin = BestMem(i_pop_number,:);
                    inite(1,i_pop_number) = 0;
                    %disp(['Best: ',num2str(fmin)])
                else
                    inite(1,i_pop_number) = inite(1,i_pop_number) + 1;
                end
                
                % What happens to the following line is BestVal is not below the fmin
                % value and therefore xmin is not defined? This does not happen - see
                % the value of fmin!!!
                xref(i_pop_number,:) = xmin;
                
                
                d_region_attraction(:,end) = [norm(ArchiveBM_LM(1,1:D,end)-ArchiveBM_LM(2,1:D,end)); 1];
                warning_LS(i_pop_number) = 0;
                
                % Check if this minimum has already been detected
                for i = 1 : size(ArchiveBM_LM,3)-1
                    
                    distance_LM_wrt_LM   = norm(ArchiveBM_LM(2,1:D,end)-ArchiveBM_LM(2,1:D,i));
                    difference_LM_wrt_LM = norm(ArchiveBM_LM(2,D+1,end)-ArchiveBM_LM(2,D+1,i));
                    
                    if (  distance_LM_wrt_LM <= norm(DELTA) * 0.001 ) && (difference_LM_wrt_LM <= norm(ArchiveBM_LM(1,D+1,end)-ArchiveBM_LM(1,D+1,i)))
                        
                        % Update dimension of the region of attraction
                        d_region_attraction(1,end) = min(d_region_attraction(1,end), d_region_attraction(1,i));
                        
                        d_region_attraction(1,i) = d_region_attraction(1,end);
                        d_region_attraction(2,i) = d_region_attraction(2,i)+1;
                        d_region_attraction(2,end) = d_region_attraction(2,end)+1;
                        
                    end
                    
                end
                
                % Define dimension of the bubble for the local restart based on
                % d_region_attraction
                % 2nd row:             % Increase the following factor by 1 every time the
                % corresponding value of exsteploc is updated (we can know how
                % precise is a bubble dimension)
                % 3rd row:             % The following factor will go to 1 when a local restart with
                % this bubble dimension will be realized
                
                
                
            end
            
        end
        
        
        % ---------------------------------------------------------------------
        % Check condition related to maximum number of function evaluations
        % ---------------------------------------------------------------------
        if sum(nFeVal) >= record(step)*nFeValMax
            for i_pop_number = 1 : pop_number
                memories(:,:,i_pop_number)    = sortrows(memories(:,:,i_pop_number),D+1);
                memories_out(step,:,i_pop_number) = memories(1,:,i_pop_number);
            end
            if sum(nFeVal)>= nFeValMax
                sum(nFeVal)
                return
            else
                step = step + 1;
            end
        end
        
        
    end
    
    
    
    % =========================================================================
    % Creation of a matrix for the adaptation of the dimension of the bubble
    % =========================================================================
    % This matrix has to be created considering the local minima found by the
    % population in two occasion:
    % - at the first local optimizer, before the first local restart
    % - after each global restart of all the populations
    % The code for now separates these two events but the aim is to merge the
    % two while and for loop that follows in a single one.
    
    while matrix_bubble_local == 0
        
        
        % DISTANCE BETWEEN POPULATIONS LOCAL MINIMA
        min_minima_distance = 10^60;
        max_minima_distance = 0;
        
        for i = 1 : pop_number-1
            
            for k = i+1 : pop_number
                
                minima_distance = norm(BestMem(i,:) - BestMem(k,:));
                min_minima_distance = min(min_minima_distance, minima_distance);
                max_minima_distance = max(max_minima_distance, minima_distance);
                
                
            end
            
        end
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
    % The follwing lines and the previous one should be appropriately joined in
    % a single group - PENSARE A COME FARLO
    
    if rem(sum(iglob), pop_number) == 0  && any(iglob)  && matrix_bubble_global == 0
        
        for i = 1 : i_pop_number
            
            % Create a matrix composed of all the known minima up to now
            minima_populations(:,:,i) = sortrows(archivebest(:,:,i),D+1);
            minima = [minima; minima_populations(1:archivebest_elements(1,i),1:D,i)];
            clear minima_populations
        end
        
        for i = 1 : size(minima,1)
            
            for k = i+1 : size(minima,1)
                
                minima_distance = norm(minima(i,:) - minima(k,:));
                min_minima_distance = min(min_minima_distance, minima_distance);
                max_minima_distance = max(max_minima_distance, minima_distance);
                
            end
            
        end
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
                p = norm(archivebest(end,1:D,i_pop_number) - archivebest(end-1,1:D,i_pop_number));
                
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
            B_mean = [B_mean, mean(B(:,1))];
        end
        
    end
    
    % BUBBLE DIMENSIONS SAMPLED FROM PARZEN DISTRIBUTION
    % Bv is a pop_number*1 matrix which, for each population elements (row)
    % contains the bubble dimension value
    Bv = parzenself_k([], B(:,1), ones(size(B(:,1))), pop_number, 'norm', 0);
    % Avoid Bv outside some boundaries?!?!
    
    
    
    % =========================================================================
    % Local and Global Restart
    % =========================================================================
    
    for i_pop_number = 1 : pop_number
        
        
        % Link each bubble dimension value to a population
        bubble_dimension(1, i_pop_number) = Bv(i_pop_number,1);
        
        if i_pop_number == 1
            bubble.pop1 = [bubble.pop1 bubble_dimension(1,i_pop_number)];
        elseif i_pop_number == 2
            bubble.pop2 = [bubble.pop2 bubble_dimension(1,i_pop_number)];
        elseif i_pop_number == 3
            bubble.pop3 = [bubble.pop3 bubble_dimension(1,i_pop_number)];
        elseif i_pop_number == 4
            bubble.pop4 = [bubble.pop4 bubble_dimension(1,i_pop_number)];
        end
        
        % Do not try to perform local or global restart if the population has
        % been already globally restarted and is waiting for the other
        % populations to be globally restarted too
        if iglob(1,i_pop_number) == min(iglob)
            
            if ~exist('warning_LS','var') || (exist('warning_LS','var') && warning_LS(i_pop_number) ~= 1)
                
                
                % --------- LOCAL RESTART
                
                % EACH population is restarted in a local bubble around xref.
                % The boundaries of the bubble are defined throught the following
                % two lines
                XVminl = max([xref(i_pop_number,:) - bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vlb]);
                XVmaxl = min([xref(i_pop_number,:) + bubble_dimension(1,i_pop_number) * DELTA/norm(DELTA); vub]);
                
                % Population re-initialized in the bubble
                pop(:,:,i_pop_number) = lhsdesign(NP,D,'criterion','maximin').*repmat(XVmaxl-XVminl,NP,1)+repmat(XVminl,NP,1);
                
                
            elseif warning_LS == 1
                
                
                % --------- GLOBAL RESTART
                
                % Count the number of global restart for each population
                iglob(1,i_pop_number) = iglob(1,i_pop_number) + 1;
                
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
                % For the definition of row_number and explanation regarding
                % the elements of archive_best considered, refer to the
                % comments below for row_number.
                %             [ClusterCenters,DataClusters,DatainClusters] = Clustering_MeanShift(archivebest(row_number+1:max(archivebest_elements),1:D,i_pop_number)',distrmin);
                [ClusterCenters,DataClusters,DatainClusters] = Clustering_MeanShift(archiveALL(:,1:D)',distrmin);
                % Number of clusters
                usc = size(ClusterCenters,2);
                
                % Center of each cluster
                meanclu=ClusterCenters';
                
                pop0=[];
                
                while size(pop0,1) < NP
                    
                    xref0 = lhsdesign(NP,D,'criterion','maximin').* repmat(vub-vlb,NP,1) + repmat(vlb,NP,1);
                    %                 xref0 = lhsdesign(NP,D,'criterion','maximin');
                    
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
                
                % The number of local restart is brought back to zero
                inite(1,i_pop_number) = 0;
                
                % PEnso che la seguente riga non serva. Vedere se e meglio
                % sostituirla per avere un codice piu elegante nella parte
                % superiore - a quanto pare no. se contraction non e` zero non
                % faccio il DE step
                contraction(1,i_pop_number) = 0;
                
                % Before starting with the next series of local restart bring
                % the fmin values back to their original values
                fmin=1e15 * ones(1,pop_number);
            end
            
        end
        % ---------------------------------------------------------------------
        
        
        
        
        %
        %    % Display information
        %    if pop_number>1 && refresh == 1
        %            if rem(sum(iglob), pop_number) ~= 0  && iglob(1,i_pop_number) > 0
        %                str_pop = int2str(i_pop_number);
        %                str = strcat('- Population',str_pop,'globally restarted and now waiting for global restart of other populations');
        %                disp(str)
        %              else
        %                str_pop = int2str(i_pop_number);
        %                str = strcat('- Local restart realized for population ',str_pop);
        %                disp(str)
        %            end
        %     end
        
        % archivebest now is a 3D matrix for each population. It collect all the
        % BestMem and BestVal values for each population, for each local or
        % global restart.
        % Since global restart are not likely to happens together for all the
        % populations, we will find ourselves in a situation in which population A
        % has already been globally restarted, but is waiting for population B to
        % be globally restarted too. While population B keeps going through local
        % optimization, it has to collect the new generated BestMem and BestVal
        % values, while population A has not such a value to collect. However,
        % since the 3D matrix does not allow to have different matrix with
        % different number of rows, we generate new rows for archivebest for
        % population A (even if they are not necessary) and fill them with NaN
        % values.
        % When both A and B will be globally restarted we will find this situation
        % for archivebest:
        % archivebest(:,:,1) = [BestMem1 BestVal1 nFeVal1;...
        %                       BestMem2 BestVal2 nFeVal2;...
        %                       BestMem3 BestVal3 nFeVal3;...
        %                       NaN      NaN      NaN;...
        %                       NaN      NaN      NaN];
        
        % archivebest(:,:,2) = [BestMem1 BestVal1 nFeVal1;...
        %                       BestMem2 BestVal2 nFeVal2;...
        %                       BestMem3 BestVal3 nFeVal3;...
        %                       BestMem4 BestVal4 nFeVal4;...
        %                       BestMem5 BestVal5 nFeVal5];
        
        % Now therefore both population have been globally restarted, and they
        % start collecting new BestMem, BestVal and nFeVal parameters in new rows
        % for each local optimization.
        % Suppose that once more population A globally restart before population B.
        % We will have the following situation for population :
        % archivebest(:,:,1) = [BestMem1 BestVal1 nFeVal1;...
        %                       BestMem2 BestVal2 nFeVal2;...
        %                       BestMem3 BestVal3 nFeVal3;...
        %                       NaN      NaN      NaN;...
        %                       NaN      NaN      NaN;...
        %                       BestMem1 BestVal1 nFeVal1;...
        %                       BestMem2 BestVal2 nFeVal2];
        % Population 1 has to perform a global restarting using cluster centers
        % that are defined using ONLY the last two rows of archivebest!!!
        % Therefore we need a number identifying the first row that we need to
        % consider for the definition of the cluster centers. This row number is
        % identified by the following parameter row_number.
        % In the previous example row_number would be 6.
        %
        % Now let us analyze how the definition of row_number is evaluated.
        % First of all we compute it only when ALL population have been globally
        % re-initialized (condition defined by the sum of the waiting vector) and
        % only if we have had a globally restart for all the population (row_number
        % is zero when the populations have not been already globally restarted at
        % least once).
        % It is computed as maximum value of archivebest_elements, which is the
        % vector which, for each population, counts how many significant (non NaN)
        % elements each population count in the archivebest matrix.
        % For the previous case it would have been
        % archivebest_elements = [3 5];
        
        if  sum(waiting) == pop_number  && iglob(i_pop_number)>0
            row_number_old = row_number;
            row_number = row_number + max(archivebest_elements);
        end
        
        
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
    
    % PLOT FOR LOCAL MINIMA AND BUBBLE DIMENSION
    if options.plot_flag == 1
        h = figure;
        S = int2str(n_figure);
        filename = S;
        aidea_movie2(archivebest,archivebest_elements,1,bubble_dimension, D, iglob, row_number_old);
        print(h, '-djpeg', filename);
        n_figure = n_figure + 1;
    end
    
    
    % We have had the first local restart, so at the next cycle we can adapt
    % the matrix for the dimension of the bubble
    first_local_restart = 0;
    
    % PLOT FOR POPULATION ADVANCEMENT - UNCOMMENT TO SEE IT
    % if options(37)==1 && i_pop_number == pop_number
    %     n_figure = aidea_movie(options,pop,n_figure);
    % end
    
    
    % ========================================================================
    % end of -------- while 1
    % that comprise both the DE step + contraction condition evaluation and the
    % local/global restart part
end
% =========================================================================

for i_pop_number = 1 : pop_number
    memories(:,:,i_pop_number) = sortrows(memories(:,:,i_pop_number),D+1);
end

% =========================================================================
% end of the function
end
% =========================================================================






%% Differential Evolution with adaptive F and CR and local search

function [BestMem, BestVal, nFeVal, pop, Val, iter, vval, new_elements] = DE_step(fname, lb, ub, pop, nFeVal, i_pop_number, mmdist, P6, archivebest, dd_limit, DE_strategy, nFeValMax, varargin)

%    INPUT
%           fname        : function handle to cost function
%           lb,ub        : lower and upper boundaries
%           pop          : population
%           nFeVal       : vector with the number of functions evaluations for
%                          each population
%           i_pop_number : number of the current population
%           mmdist       : convergence threshold (for population contraction)
%           P6           : probability of running selected DE strategies
%           archivebest  :
%           dd_limit     :
%           DE_strategy  : integer to identify selected DE strategies
%           nFeValMax    : maximum number of function evaluations allowed
%           varargin     : additional inputs
%
%   OUTPUT
%           BestMem       : best solution vector
%           BestVal       : best solution value
%           nFeVal        : number of function evaluations
%           pop           : population
%           Val           : cost function associated to pop
%           iter          :
%           vval          :
%           new_elements  :



% (c) Edmondo Minisci and Massimiliano Vasile 2013
%     Marilena Di Carlo, 2015

% =========================================================================
% Initialization
% =========================================================================
% ?
new_elements = 0;

% ?
step_DE = 0;

% The population is composed by NP elements with dimension D
[NP,D]    = size(pop);

% Val is what in idea2.m was called fitness
Val       = zeros(1,NP);

% Best population member ever
BestMem   = zeros(1,D);

% Mean element of the population
MeanPop   = mean(pop);

% Maximum distance between elements of the population and their mean value
mmdistm = -1.e10;

for imm = 1 : size(pop)
    dista = norm(pop(imm,:)-MeanPop);
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


% Number of populations
pop_number = size(archivebest,3);

% =========================================================================
% Function Evaluation for Initial Population (Parents)
% =========================================================================
% start with first population member
ibest   = 1;

Val(1)  = feval(fname,pop(ibest,:)',varargin{:});
BestVal = Val(1);                 % best objective function value so far
nFeVal(1,i_pop_number)  = nFeVal(1,i_pop_number) + 1;

if sum(nFeVal) >= nFeValMax
    vval = 0;
    iter = 0;
    new_elements = 1;
    return
end

for i = 2 : NP                        % check the remaining members
    Val(i) = feval(fname,pop(i,:)',varargin{:});
    nFeVal(1,i_pop_number) = nFeVal(1,i_pop_number) + 1;
    
    if (Val(i) < BestVal)           % if member is better
        ibest   = i;                 % save its location
        BestVal = Val(i);
    end
    
    if sum(nFeVal) >= nFeValMax
        new_elements = i;
        vval = 0;
        iter = 0;
        return
    end
    
end

CurrBest = pop(ibest,:);            % best member of current iteration
BestMem = CurrBest;                 % best member ever



% =========================================================================
% DE Initialization
% =========================================================================

BestMat  = zeros(NP,D);                 % initialize bestmember  matrix
InterPop  = zeros(NP,D);                % intermediate population of perturbed vectors

RotIndArr = (0:1:NP-1);               % rotating index array (size NP)

RotIndArr2  = zeros(NP);                % another rotating index array

IndArr1  = zeros(NP);                % index arrays
IndArr2  = zeros(NP);
IndArr3  = zeros(NP);
IndArr4  = zeros(NP);
IndArr5  = zeros(NP);
IndPointArr = zeros(4);

iter = 1;

mmdistm=1.e10;

nostop = 1;
% nFeValmem = nFeVal;
% itermem = 0;


%% Initialization of CRF to uniform distribution

% Crossover probability CR is defined in the interval 0.1 to 0.99
% CRa is a (D+1)*3 matrix
% If D = 3 CRa will be
% CRa = [0.1  0  0;...
%        CR2  0  0;...
%        CR3  0  0;...
%        0.99 0  0];

% delta     = (.99-.1)/(D);
delta = (.99-.1)/(D/2);
CRa_1st   = (0.1:delta:0.99)';          % 1st column of matrix CRa
CRa       = [CRa_1st zeros(size(CRa_1st)) zeros(size(CRa_1st))];

% Differential weight F is defined in the interval -1 to 1
% Fa is a (D+1)*3 matrix
% delta     = (1-(-1))/(D);
delta     = (1-(-1))/(D/2);
Fa_1st    = (-1:delta:1)';             % 1st column of matrix Fa
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
% and therefore its dimensions will be (D+1)(D+1)x4
for im = 1 : size(CRa(:,1))
    for in = 1 : size(Fa(:,1))
        CRFa=[CRFa; CRa(im,1) Fa(in,1) 0 0];
    end
end

while nostop
    
    % =====================================================================
    % CR and F are sampled from the Parzen distribution
    % =====================================================================
    % CRFv is a NP*2 matrix which, for each population elements (row)
    % contains the CR (1st column) and F (2nd column) values
    CRFv=parzenself_k([],CRFa(:,1:2),ones(size(CRFa(:,1))),NP,'norm',0);
    
    % Separate CR from F values
    % At the end of the following lines:
    % CRv -> pop_number*1 matrix
    % Fv  -> pop_number*1 matrix
    CRv=CRFv(:,1);
    Fv=CRFv(:,2);
    
    % Avoid CR outside the boundaries 0.1 to 0.99
    % If some of the obtained CRv associated to a population element is
    % lower than 0.1, its values is brought back to 0.1
    % At the same time, is some element have a value bigger than 0.99, it
    % is brought back to 0.99
    % The matrix CRv as defined in the following lines is indeed
    % CRv = [CR1 CR2 CR3 ..... CRn;...
    %        0.1 0.1 0.1       0.1]
    % where CR1, CR2 etc represent the CR values associated to the 1st, 2nd
    % and so on element of the population
    % The max function returns a row vector containing the maximum value in
    % each column
    % The same holds for the min comparison
    % At the end of the following two lines CRv will be a pop_numberx1 column vector
    % (it is transposed back at the end of the lines!)
    CRv=(max([CRv'; 0.1*ones(size(CRv'))]))';
    CRv=(min([CRv';  0.99*ones(size(CRv'))]))';
    CR=repmat(CRv,1,D);
    
    Fv=(max([Fv'; -1*ones(size(Fv'))]))';
    Fv=(min([Fv';  1*ones(size(Fv'))]))';
    F=repmat(Fv,1,D);
    
    
    
    % =========================================================================
    % Offspring generation
    % =========================================================================
    
    % =====================================================================
    % Probabilistic selection of the strategy (and count how many times a
    % given strategy is applied)
    % =====================================================================
    if rand(1,1) < P6
        %
        strategy = 1;
    else
        
        strategy = 2;
    end
    
    
    
    % =====================================================================
    % Creation of the offspring population
    % =====================================================================
    
    popold = pop;                   % save the old population
    
    IndPointArr = randperm(4);              % index pointer array
    
    IndArr1  = randperm(NP);             % shuffle locations of vectors
    RotIndArr2 = rem(RotIndArr+IndPointArr(1),NP);        % rotate indices by IndPointArr(1) positions
    IndArr2  = IndArr1(RotIndArr2+1);                 % rotate vector locations
    RotIndArr2 = rem(RotIndArr+IndPointArr(2),NP);
    IndArr3  = IndArr2(RotIndArr2+1);
    RotIndArr2 = rem(RotIndArr+IndPointArr(3),NP);
    IndArr4  = IndArr3(RotIndArr2+1);
    RotIndArr2 = rem(RotIndArr+IndPointArr(4),NP);
    IndArr5  = IndArr4(RotIndArr2+1);
    
    PopMat1 = popold(IndArr1,:);             % shuffled population 1
    PopMat2 = popold(IndArr2,:);             % shuffled population 2
    PopMat3 = popold(IndArr3,:);             % shuffled population 3
    PopMat4 = popold(IndArr4,:);             % shuffled population 4
    PopMat5 = popold(IndArr5,:);             % shuffled population 5
    
    for i=1:NP                      % population filled with the best member
        BestMat(i,:) = CurrBest;          % of the last iteration
    end
    
    
    st = strategy;		                     % binomial crossover
    % Mask e
    % To each element of the population is associated a different value of
    % CR!
    MaskInterPop = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise
    MaskOldPop = MaskInterPop < 0.5;         % inverse mask to MaskInterPop
    
    
    % According to the randomly selected value of strategy, the DE will
    % produce a mutant vector using a completely random vector or a random
    % vector taken from the archive of minima of all the populations or
    % from the best member of the population
    
    % ---------------------------------------------------------------------
    % DE_strategy = 1 : DE/best, DE/rand
    % ---------------------------------------------------------------------
    
    if DE_strategy == 1
        
        if (st == 1)                                                            % DE/best/1
            InterPop = popold  + F.*(BestMat - popold) + F.*(PopMat1 - PopMat2);
            %             InterPop = BestMat + F.*(PopMat1 - PopMat2);                        % differential variation
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;             % crossover
        elseif (st == 2)                                                        % DE/rand/1
            InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
        end
        
        % ---------------------------------------------------------------------
        % DE_strategy = 3 : DE/best, DE/arch
        % ---------------------------------------------------------------------
    elseif DE_strategy == 3
        
        if (st == 1)                                                            % DE/best/1
            InterPop = BestMat + F.*(PopMat1 - PopMat2);                        % differential variation
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;             % crossover
        elseif (st == 2)                                                        % DE/rand/1
            % When the population have not yet reached condition of contraction
            % archivebest is still empty and therefore a random vector is
            % chosen instead. Otherwise use DE/best here!!
            if isempty(archivebest)
                InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            else
                
                while 1
                    
                    % Define a random index to choose which archive population
                    % to use
                    archive_number = ceil(rand*pop_number);
                    
                    % Define a random index for the element of archivebest
                    random_index = ceil(rand*size(archivebest,1));
                    
                    % archivebest is composed of NaN values (see above). We can
                    % not accept these values therefore we perform a check on
                    % the selected value of archivebest. If the selected values
                    % is acceptable we assign it to PopMat3
                    if isnan(archivebest(random_index,1,archive_number))
                    else
                        for i = 1 : NP
                            PopMat3(i,1:D) = archivebest(random_index,1:D,archive_number);
                        end
                        break
                    end
                end
                InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            end
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
            
        end
        
        % ---------------------------------------------------------------------
        % DE_strategy = 2 : DE/rand, DE/arch
        % ---------------------------------------------------------------------
    elseif DE_strategy == 2
        
        if (st == 1)                                                            % DE/best/1
            % When the population have not yet reached condition of contraction
            % archivebest is still empty and therefore a random vector is
            % chosen instead. Otherwise use DE/best here!!
            if isempty(archivebest)
                InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            else
                
                while 1
                    
                    % Define a random index to choose which archive population
                    % to use
                    archive_number = ceil(rand*pop_number);
                    
                    % Define a random index for the element of archivebest
                    random_index = ceil(rand*size(archivebest,1));
                    
                    % archivebest is composed of NaN values (see above). We can
                    % not accept these values therefore we perform a check on
                    % the selected value of archivebest. If the selected values
                    % is acceptable we assign it to PopMat3
                    if isnan(archivebest(random_index,1,archive_number))
                    else
                        for i = 1 : NP
                            PopMat3(i,1:D) = archivebest(random_index,1:D,archive_number);
                        end
                        break
                    end
                end
                InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            end
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
            
        elseif (st == 2)                                                        % DE/rand/1
            InterPop = PopMat3 + F.*(PopMat1 - PopMat2);
            InterPop = popold.*MaskOldPop + InterPop.*MaskInterPop;
        end
        
    end
    
    % =====================================================================
    % Select which vectors are allowed to enter the new population
    %======================================================================
    CRFa=sortrows(CRFa, [3 4]);
    
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
        % Every component violating the boundaries is projected back into D
        % by picking a new value dependent upon the boundaries and a random
        % uniform distribution (pag.5)
        % =================================================================
        
        % Distance of the i-th element of the population from the lower
        % boundaries
        dInterPopl = InterPop(i,:)-lb;
        
        % How many component of the i-th element of the population are
        % below the lower boundary?
        ElTooLow  = find(dInterPopl<0);
        lElTooLow = length(ElTooLow);
        
        
        % Distance of the i-th element of the population from the upper
        % boundaries
        dInterPopr = ub-InterPop(i,:);
        
        % How many component of the i-th element of the population are
        % above the lower boundary?
        ElTooHigh  = find(dInterPopr<0);
        lElTooHigh = length(ElTooHigh);
        
        
        % Bring back the considered element into the boundaries of the
        % problem
        if lElTooLow > 0
            %             InterPop(i,ElTooLow) = lb(ElTooLow) + rand(1,lElTooLow).*(ub(ElTooLow) - lb(ElTooLow));
            InterPop(i,ElTooLow) = ( popold(ElTooLow) + lb(ElTooLow) ) / 2;
        end
        if lElTooHigh > 0
            %             InterPop(i,ElTooHigh) = lb(ElTooHigh) + rand(1,lElTooHigh).*(ub(ElTooHigh) - lb(ElTooHigh));
            InterPop(i,ElTooHigh) = ( popold(ElTooHigh) + ub(ElTooHigh) ) / 2;
        end
        
        
        % =================================================================
        % Value of the function f for each new element of the population
        % =================================================================
        
        TempVal = feval(fname,InterPop(i,:)',varargin{:});   % check cost of competitor
        nFeVal(1,i_pop_number)  = nFeVal(1,i_pop_number) + 1;
        
        
        if (TempVal <= Val(i))  % if competitor is better than value in "cost array"
            
            
            dd=abs((TempVal-Val(i)));
            
            
            % replace old vector with new one (for new iteration)
            pop(i,:) = InterPop(i,:);
            
            % save value in "cost array"
            Val(i) = TempVal;
            
            
            % Adattamento di CR e F nella popolazione
            for ic=1:size(CRFa,1)
                
                % If the previous difference between parent and offspring
                % function value is smaller than the new one (under the
                % condition that the offspring perform better than the
                % parent - we are already in an if condition, that is if
                % TempVal <= Val(i) )
                if CRFa(ic,3)<dd
                    
                    % The substitution of CR is subject to another
                    % limitation
                    if abs(dd)>dd_limit
                        
                        % Substitution of CR
                        CRFa(ic,1)=CRv(i);
                    end
                    
                    % Substitue the previous value of F with the one used
                    % for the considered offspring
                    CRFa(ic,2)=Fv(i);
                    
                    % Associate the corresponding dd and iter value
                    CRFa(ic,3)=dd;
                    CRFa(ic,4)=iter;
                    break
                end
            end
            
            
            if (TempVal < BestVal)     % if competitor better than the best one ever
                BestVal = TempVal;      % new best value
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
        dista = norm(pop(imm,:)-MeanPop);
        if dista > mmdistm
            mmdistm=dista;
        end
    end
    
    % mmdistm is the maximum distance between population mean and elements
    % of the population
    %[mmdistm min(val)]
    %
    %     if isave < size(vsave,2)
    %         if nFeVal >= vsave(1,isave+1)
    %             isave = isave + 1;
    %             savesol = [savesol; BestMem BestVal nFeVal];
    %         end
    %     end
    %
    % Collects results of the considered offspring
    vvalr = [iter min(Val) mean(Val) max(Val) mmdistm mmdistm/mmdistm0];
    vval = [vval; vvalr];
    %
    iter = iter + 1;
    %
    %     % nFeVal is equal to nFeValmem only if the population has not advanced
    %     % from parent to offspring. If the population has not advanced, itermem
    %     % increases by one
    %     if nFeVal~=nFeValmem
    %         nFeValmem=nFeVal;
    %         itermem=0;
    %     else
    %         itermem=itermem+1;
    %     end
    
    step_DE = step_DE + 1;
    
    nostop=mmdistm/max(vval(:,5))>mmdist  && sum(nFeVal)<nFeValMax && step_DE < (D/10)*100 ;
    % nostop=mmdistm/max(vval(:,5))>mmdist  && sum(nFeVal)<nFeValMax  ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     disp(['DE generation: ' num2str(step_DE)]);
    % %     disp(strcat('\rho for DE convergence: ', num2str(mmdistm/max(vval(:,5)) - mmdist)))
    %     disp(['Minimum: ' num2str(BestVal)]);
    %     disp('------------------------------------------------')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end



