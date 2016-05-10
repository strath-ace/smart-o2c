function [memory,nfeval]=macs7v15c2_old(func,memory,vlb,vub,options,filename,fileload,varargin)

%  memory=macs7v15c2(func,vlb,vub,options,filename,fileload,varargin)
%
%  INPUT
%         func    : function handle
%         vlb,vub : boundaries of the search domain, as vectors
%         options : structure of optimization parameters
%                   'maxnfeval' max num of function evaluations albeit
%                               forced in a "weak" way (default 5000)
%                   'popsize'   max num of agents (default 10)
%                   'rhoini'    Initial local hypercube size (default 1)
%                   'F'         F parameter for Differential Evolution
%                               (default 0.9)
%                   'CR'        crossover probability (default 0.9)
%                   'p_social'  Ratio between population performing social 
%                               moves and total population (default 0.2)
%                   'max_arch'  max size of the archive at output (default
%                               10)
%                   'coord_ratio'ratio of coordinates to check after which 
%                               no further search will be performed, as in 
%                               Zuiani 3b, pag 6 (default 0.25)
%                   'contr_ratio' Contraction ratio for rho (default 0.5)
%                   'draw_flag' plot flag (default 0)
%                   'cp'        constraint flag (default 0)
%                   'MBHflag'   MBH steps flag (default 0)
%                               Currently only for unconstrained
%                   'cpat'      pattern to DE vs pattern to local flag
%                               (default 0)
%                   'DE_strategy' Strategy to use for DE (ie pull towards
%                               best element or random one, default 'best')
%                   'explore_all' flag to state if all agents must perform
%                               local exploration or not (default 0)
%
%  OUTPUT
%           memory :  matrix containing the archived solutions and the
%                       associated value of the cost function
%                       memory=[x f]
%           nfeval   :  number of total function evaluations

% CHANGE LOG
% Created by: Federico Zuiani 2011
% Revised, cleaned and optimized by: Lorenzo A. Ricciardi 2015

%%  MACS PARAMETERS DEFAULT VALUES

default.maxnfeval=5000;
default.popsize=10;
default.min_popsize=3;
default.rhoini=1;
default.F=0.9;
default.Fmax=2;             %upper bound to F
default.CR=0.9;
default.p_social=0.2;
default.max_arch=10;
default.coord_ratio=0.25;
default.contr_ratio=0.5;
default.draw_flag=0;
default.cp=0;
default.MBHflag=0;
default.cpat=0;
default.explore_DE_strategy='best';
default.social_DE_strategy='DE/current-to-rand/1';
default.explore_all=1;
default.v=0;
default.int_arch_mult=1.5;

%%  MACS PARAMETERS VALIDATION AND DEFAULTING

if isfield(options,'maxnfeval')
   if ~isPositiveIntegerValuedNumeric(options.maxnfeval)
       error('maxnfeval must be a positive integer!');
   end
else
    warning(['maxnfeval not supplied, using default value ' num2str(default.maxnfeval)]);
    options.maxnfeval = default.maxnfeval;
end

if isfield(options,'popsize')
   if ~isPositiveIntegerValuedNumeric(options.popsize)
       error('popsize must be a positive integer!');
   else
       if (options.popsize<default.min_popsize)
            error('popsize must be at least 3!');
       end
   end
else
    warning(['popsize not supplied, using default value ' num2str(default.popsize)]);
    options.popsize = default.popsize;
end

if isfield(options,'rhoini')
   if ~isnumeric(options.rhoini)
       error('rhoini must be a real number between 0 and 1')
   else
       if (options.rhoini<0)||(options.rhoini>1)
           error('rhoini must be between 0 and 1!');
       end
   end
else
    warning(['rhoini not supplied, using default value ' num2str(default.rhoini)]);
    options.rhoini = default.rhoini;
end

if isfield(options,'F')
   if ~isnumeric(options.F)
       error('F must be a real number between 0 and 2')
   else
       if (options.F<0)||(options.F>default.Fmax)
           error(['F must be between 0 and ' num2str(default.Fmax) '!']);
       end
   end
else
    warning(['F not supplied, using default value ' num2str(default.F)])
    options.F = default.F;
end

if isfield(options,'CR')
   if ~isnumeric(options.CR)
       error('CR must be a real number between 0 and 1')
   else
       if (options.CR<0)||(options.CR>1)
           error('CR must be between 0 and 1!');
       end
   end
else
    warning(['CR not supplied, using default value ' num2str(default.CR)])
    options.CR = default.CR;
end

if isfield(options,'p_social')
   if ~isnumeric(options.p_social)
       error('p_social must be a real number between 0 and 1')
   else
       if (options.p_social<0)||(options.p_social>1)
           error('p_social must be between 0 and 1!');
       end
   end
else
    warning(['p_social not supplied, using default value ' num2str(default.p_social)]);
    options.p_social = default.p_social;
end

if isfield(options,'max_arch')
   if ~isPositiveIntegerValuedNumeric(options.max_arch)
       error('max_arch must be a positive integer!');
   end
else
    warning(['max_arch not supplied, using default value ' num2str(default.max_arch)]);
    options.max_arch = default.max_arch;
end

if isfield(options,'coord_ratio')
   if ~isnumeric(options.coord_ratio)
       error('coord_ratio must be a real number between 0 and 1')
   else
       if (options.coord_ratio<0)||(options.coord_ratio>1)
           error('coord_ratio must be between 0 and 1!');
       end
   end
else
    warning(['coord_ratio not supplied, using default value ' num2str(default.coord_ratio)]);
    options.coord_ratio = default.coord_ratio;
end

if isfield(options,'contr_ratio')
   if ~isnumeric(options.contr_ratio)
       error('contr_ratio must be a real number between 0 and 1')
   else
       if (options.contr_ratio<0)||(options.contr_ratio>1)
           error('contr_ratio must be between 0 and 1!');
       end
   end
else
    warning(['contr_ratio not supplied, using default value ' num2str(default.contr_ratio)]);
    options.contr_ratio = default.contr_ratio;
end

if isfield(options,'draw_flag')
   if ~isPositiveIntegerValuedNumeric(options.draw_flag+1) %0 is not counted as postive, this should do the trick
       error('draw_flag must be a positive integer!');
   end
else
    warning(['draw_flag not supplied, using default value ' num2str(default.draw_flag)]);
    options.draw_flag = default.draw_flag;
end

if isfield(options,'cp')
   if ~isPositiveIntegerValuedNumeric(options.cp+1) %0 is not counted as postive, this should do the trick
       error('cp must be a 0 or 1!');
   else
       if options.cp>1
           error('cp must be a 0 or 1!');           
       end
   end
else
    warning(['cp not supplied, using default value ' num2str(default.cp)]);
    options.cp = default.cp;
end

if isfield(options,'MBHflag')
   if ~isPositiveIntegerValuedNumeric(options.MBHflag+1) %0 is not counted as postive, this should do the trick
       error('MBHflag must be a 0 or 1!');
   else
       if options.MBHflag>1
           error('MBHflag must be a 0 or 1!');           
       end
   end
else
    warning(['MBHflag not supplied, using default value ' num2str(default.MBHflag)]);
    options.MBHflag = default.MBHflag;
end

if isfield(options,'cpat')
   if ~isPositiveIntegerValuedNumeric(options.cpat+1) %0 is not counted as postive, this should do the trick
       error('cpat must be a 0 or 1!');
   else
       if options.cpat>1
           error('cpat must be a 0 or 1!');           
       end
   end
else
    warning(['cpat not supplied, using default value ' num2str(default.cpat)]);
    options.cpat = default.cpat;
end

if isfield(options,'explore_DE_strategy')
   if ~(strcmp(options.explore_DE_strategy,'best') || strcmp(options.explore_DE_strategy,'rand'))
      error('Unrecognized explore_DE_strategy: must be either best or rand ') 
   end
else
    warning(['explore_DE_strategy not supplied, using defauld value ' default.explore_DE_strategy])
    options.explore_DE_strategy = default.explore_DE_strategy;
end

if isfield(options,'social_DE_strategy')
   if ~(strcmp(options.social_DE_strategy,'DE/rand/1') || strcmp(options.social_DE_strategy,'DE/current-to-rand/1'))
      error('Unrecognized social_DE_strategy: must be either DE/rand/1 or DE/current-to-rand/1') 
   end
else
    warning(['social_DE_strategy not supplied, using defauld value ' default.social_DE_strategy])
    options.social_DE_strategy = default.social_DE_strategy;
end

if isfield(options,'explore_all')
   if ~isPositiveIntegerValuedNumeric(options.explore_all+1) %0 is not counted as postive, this should do the trick
       error('explore_all must be a 0 or 1!');
   else
       if options.explore_all>1
           error('explore_all must be a 0 or 1!');           
       end
   end    
else
   warning(['explore_all not supplied, using default value ' default.explore_all])
   options.explore_all = default.explore_all;
end

if isfield(options,'v')
   if ~isPositiveIntegerValuedNumeric(options.v+1) %0 is not counted as postive, this should do the trick
       error('v must be a 0 or 1!');
   else
       if options.v>1
           error('v must be a 0 or 1!');           
       end
   end    
else
   warning(['v not supplied, using default value ' default.v])
   options.v = default.v;
end

if isfield(options,'int_arch_mult')
   if ~isnumeric(options.int_arch_mult)
       error('int_arch_mult must be a real number greater than 1')
   else
       if (options.int_arch_mult<1)
           error('int_arch_mult must be a real number greater than 1');
       end
   end
else
    warning(['int_arch_mult not supplied, using default value ' num2str(default.int_arch_mult)]);
    options.int_arch_mult = default.int_arch_mult;
end

%%  MACS SUBPARAMETERS AUTOSETTING

n_social        = round(options.p_social*options.popsize);                  % number of agents performing social actions
supbr_upd_freq = n_social;                                                  % number of iterations after which new subproblems are chosen
lx             = length(vlb);                                               % number of parameters of objective functions

%%  MEMORY INITIALIZATION, ONLY IF MEMORY IS EMPTY

if isempty(memory)                                                          % if memory is already populated, no initialization is needed!
    
    %% CHOOSE RANDOM PARAMETERS AND SET UP PROBLEM DIMENSIONALITY AND CONSTRAINTS
    
    x=lhsu(vlb,vub,options.popsize);                                        % latin hypercube sampling over the whole domain
    maxC=zeros(options.popsize,1);                                          % maximum constraint violation
    
    % Check number of objectives and constraints
    if options.cp==0                                                        % if problem is unconstrained
        
        [f_dummy]=func(x(1,:),varargin{:});                                 % make a bogus evaluation of the fcn
        mfit=length(f_dummy);                                               % just to measure it's dimensions
        ncon=1;                                                             % and set number of constraints to 1, will be used later
        
    else                                                                    % if proble is constrained
        
        [f_dummy,c_dummy]=func(x(1,:),varargin{:});                         % do the same as above
        mfit=length(f_dummy);
        ncon=length(c_dummy);                                               % but set the correct number of constraints, to be used later
        
    end
    
    first_shrink=1;                                                         % shrink archive just after initialization, archive will contain best first candidate solutions
    c = zeros(options.popsize,ncon);                                        % initialize the matrix of constraint violations
    
    %% Prototype adaptiviy for CR -> use the utility function to adapt in the main loop (NOT USED, YET)
    %           CR=ones(1,2);
    %                 idCR=dprob(CR/max(CR));
    %                 c=idCR/100;
    %              CR(idCR)=CR(idCR)+2*(impr-0.5);
    %            CR=CR-min(CR);
    %
    
    %% INITIALIZATION OF LAMBDA VECTORS DEFINING SUBPROBLEMS
    
    switch mfit
        
        case 1                                                              % Single Objective
            
            n_lambda=2;                                                     % number of subproblems
            tmp_max_arch=round(options.max_arch*options.int_arch_mult);      % temporary archive size, used to store the population just after initialization
            lambda=[1;1];                                                   % at least 2 subproblems are needed for the whole strategy to work
            
        case 2                                                              % Bi-objective
            
            n_lambda=mfit*100;                                              % number of subproblems is 100 times the number of objectives, so always 200... that 100 could be a parameter...
            tmp_max_arch=round(min([n_lambda options.max_arch])*options.int_arch_mult);       % number of subproblems is the greater of 1.5 times the num of subproblems or 1.5 times the max size of the archive
            alfas=linspace(0,pi/2,n_lambda)';                               % subproblems of a Bi-objective optimization are distributed regularly around a 1/4 circle, so here we compute the angles
            alfas=alfas(2:end-1);                                           % orthogonal subproblems (alpha=0 and pi/2) are excluded from this list, since they will be introduced manually later
            alfas=alfas(randperm(n_lambda-2));                              % shuffle the alphas, to shuffle the subproblems (why is this needed???)
            lambda=[eye(mfit); sin(alfas) cos(alfas)];                      % finally, in lambdas are m orthogonal subproblems, plus all the other ones around a 1/4 circle
            
        case 3                                                              % Tri-objective
            
            n_lambda=mfit*100;                                              % number of subproblems is 100 times the number of objectives, so always 300... that 100 could be a parameter...
            tmp_max_arch=round(min([n_lambda options.max_arch])*options.int_arch_mult);       % number of subproblems is the greater of 1.5 times the num of subproblems or 1.5 times the max size of the archive
            n_alpha=round(sqrt(n_lambda));                                  % subproblems of a Tri-objective optimization are distributed regularly around a 1/8 sphere, so here we compute the angles (always 17)
            n_beta=round(n_lambda/n_alpha);                                 % always 18
            alphas=linspace(0,1/4,n_alpha)'*2*pi;                           % compute distribution of alphas
            betas=acos(linspace(-0.5,0.5,n_beta)'-0.5)-pi/2;                % betas
            lambda=eye(3);                                                  % orthogonal subproblems
            for i=1:n_alpha                                                 % compute points on the sphere
                for j=1:n_beta
                    if (i*j~=1)&&~(i==n_alpha&&j==1)&&(i*j~=n_alpha*n_beta) % excluding the ones on the borders...
                        
                        lambda=[lambda;cos(alphas(i))*cos(betas(j)),sin(alphas(i))*cos(betas(j)),sin(betas(j))];
                    end
                end
            end
            n_lambda=n_alpha*n_beta;                                        % evaluates to 306, not 300...
            lambda(4:end,:)=lambda(3+randperm(n_lambda-3),:);
            
        otherwise                                                           % Many objectives
            
            n_lambda=mfit*100;                                              % latin hypercube sampling
            tmp_max_arch=round(min([n_lambda options.max_arch])*options.int_arch_mult);
            lambda=[eye(mfit); lhsu(zeros(1,mfit),ones(1,mfit),n_lambda-mfit)];
    end
    
    %% NORMALZATION OF VECTORS ASSOCIATED TO THE SUBPROBLEMS
    
    lambda(lambda(mfit+1:end,:)==0)=0.0001;    
    lambda=lambda./repmat(sum(lambda.*lambda,2).^0.5,1,mfit);    
    
    %% INITIALIZATION OF PIGR, V AND RHO
    
    pigr=ones(1,n_lambda);              % utility function
    v=zeros(options.popsize,lx);                % PROBABLY THIS IS THE VECTOR OF IMPROVEMENT SINCE THE LAST ITERATION    
    rho(:,1)=options.rhoini*ones(options.popsize,1);    %the radius of the local search
    rho(:,2)=zeros(options.popsize,1);          %counter of iterations, after some will reset radius to initial value
    
    %% UNUSED
    
    %trick to evaluate an unknown function with the right amount of parameters
    %6 is the number of chars before the varargin variable
    %if nargin>6&&~isempty(varargin)
    %    fevalstr='func(y,varargin{:})';
    %else
    %    fevalstr='func(y)';
    %end
    
    %% CREATION OF INITIAL POPULATION 
    
    nfeval=0;
    
    if nargin>5&&exist(fileload,'file')                                     % if there are previous evaluations of the fcn and they are stored in a file, open it
        
        load(fileload,'memory','x','f')                                   % if values are stored in the file and these values are more than the size of the population
        
        if exist('f','var')&&(size(f,1)>=options.popsize)                   % if there are more function evaluations than the max needed (popsize)
        
            x=x(1:options.popsize,:);                                       % store only popsize values in x
            f=f(1:options.popsize,:);                                       % and in f
            
        else                                                                % if not (BEWARE, they should EXIST BUT BE LESS OR EQUAL THAN POPSIZE)
            
            ids=randperm(size(memory,1));                                   % shuffle the permutation of the components of memory
            ids=ids(1:options.popsize);                                     % take only popsize of them
            x=memory(ids,1:lx);                                             % and put those items as x
            f=memory(ids,lx+1:lx+mfit);                                     % and f
        end
        
        n_mem=size(memory,1);                                               % update current number of memory in the archive
        nfeval=n_mem;                                                       % and keep in mind that the function has already been evaluated n_mem times, in the past

    else                                                                    % if no file exists
        
        f=zeros(options.popsize,mfit);                                      % initialize a matrix of zeros as f
                                        
        for i=1:options.popsize                                             % and populate it:
            
            y=x(i,:);                                                       % create a temporary y
            
            if options.cp==0                                                % if the problem is unconstrained
            
                f(i,:)=func(y,varargin{:});                                 % evaluate f
                maxC(i)=0;                                                  % and say all the constraints are satisfied
            
            else                                                            % if the problem is constrained
                
                [f(i,:),c(i,:)]=func(y,varargin{:});                        % evaluate f
                maxC(i)=max(c(i,:));                                        % and say maximum violation is max of c
            
            end
            
            nfeval=nfeval+1;                                                % then update nfeval
            
        end
        
        dom=dominance(f,0);                                                 % once you have f, evaluate dominances. BEWARE, THIS SHOULD ALSO BE DONE FOR LOADED PREVIOUS EVALUATIONS!!!
        cid=maxC;                                                           % NOT SURE IF THIS IS RIGHT... IS IS SUPPOSED TO BE THE ID OF THE MAX CONSTRAINT VIOLATION?
        
        j=0;
        
        for i=1:options.popsize                                             % endow memory with dominance and max constraint violation
            %if cid(i)<=0                                                   SEEMS THAT CONSTRAINTS WHERE ONLY SKETCHED...
            j=j+1;
            memory(j,:)=[ x(i,:) f(i,:) dom(i) cid(i) ];
            %end
        end
        %n_mem=size(memory,1);     %update memory size, number of
    end
    
    delta=max(memory(:,lx+1:lx+mfit),[],1)-min(memory(:,lx+1:lx+mfit),[],1);    % compute delta, the excursion in criteria space for each element of the population(can also be 0!)
    z=min(memory(:,lx+1:lx+mfit),[],1);                                   % and z, the min in criteria space NOTE, THIS IS NOT THE TRUE MIN, ONLY THE BEST COMPUTED SO FAR!!!!
    
    archsize=length(memory(:,1));  
        
    %% CHECK WHICH IS THE BEST APPROXIMATION OF THE I-TH SUBPROBLEM
    
    lambda_f_best=zeros(n_lambda,mfit);                                     % initialize lambda_f (same as criteria space dimensionality)
    lambda_x_best=zeros(n_lambda,lx);                                       % initialize lambda_x (same as parameter space dimensionality)

    for i=1:n_lambda                                                        % for each subproblem
    
        g_tmp = hp_g_fun(memory(:,lx+1:lx+mfit),lambda(i,:),z);             % evaluate g_fun on all agents for current i-th subproblem        
        [~,tmp_id]=min(g_tmp);                                              % store the minimum g and it's id (position), so tmp_id is the ID of the agent that best solves i-th subproblem
        lambda_f_best(i,:)=memory(tmp_id,lx+1:lx+mfit);                     % populate lambda_f_best with the f value of tmp agent
        lambda_x_best(i,:)=memory(tmp_id,1:lx);                             % populate lambda_x_best with the x value of tmp agent
    
    end
        
    %% SELECT THE INITIAL POPULATION OF SUBPROBLEMS TO BE TREATED
    
    % This complex section chooses the agents which are in a good position
    % for a subproblem, and forces them to perform SOCIAL actions with the
    % aim of improving those subproblems.
    
    if n_social>=mfit                                                       % if the number of agent performing SOCIAL actions is equal or greater than the number of objectives (so we can explore all functions at least individually)
        
        if mfit==1                                                          % if single objective
            
            act_subpr=1:n_lambda;                                           % active problems IDs are all problems (which are 2 times the same problem, for the code to work)
            
        else                                                                % if multiobjective
            
            act_subpr=randperm(n_lambda-mfit)+mfit;                         % create a vector with the (shuffled) IDs of subproblems to solve, excluding the first orthogonal ones
            act_subpr=[1:mfit act_subpr(1:n_social-mfit)];                  % the vector of active subproblems is then composed, first, of the orthogonal subproblems, and then, by an amount of random non-orthogonal subproblems equal to the remaining number of agents (i.e that will not be devoted to solve orthogonal subproblems)
                                                                            % THUS WE ALWAYS INITIALIZE N_SOCIAL SUBPROBLEMS, ONE FOR EVERY SOCIAL AGENT!!!                                                                            
        end
        
        if archsize>1                                                       % if there are at least 2 items in the initial archive (should be ALWAYS)
        
            id_pop_act_subpr=zeros(1,n_social);                             % initialize a vector which will contain the IDs of the agents which will be performing SOCIAL actions
            id_local=1:options.popsize;                                     % and a vector with the ids of the elements which will be performing LOCAL actions (initially all)
                                    
            for i=1:n_social                                                % for every agent performing SOCIAL actions
            
                tmp=Inf;                                                    % ugly init to Inf this temp variable
                
                for j=id_local                                             % for every ID of agents performing SOCIAL actions (BEWARE, THIS IS A LIST OF IDS NOT AN ORDERED LIST FROM 1 TO N_SOCIAL)
                
                    if isnan(norm((lambda_f_best(i,:)-f(j,:))./delta))      % check lambda_f_best or f contain NaN, or if they are equal and are divided by zero
                    
                        error('The objective vector contains NaN or has no variation')  % in that case, nothing can be done
                    
                    end
                    
                    if norm((lambda_f_best(i,:)-f(j,:))./delta)<tmp         % if the relative difference has diminished since last couple checked (starting from Inf at first, so sure for first iteration)
                    
                        tmp=norm((lambda_f_best(i,:)-f(j,:))./delta);       % set this diminishing as the new threshold for the next check
                        id_pop_act_subpr(i)=j;                              % and mark the ID of this (j-th) SOCIAL agent as the CLOSEST to the LOCAL (i-th) agent
                        
                    end
                    
                end          
                
                x(id_pop_act_subpr(i),:)=lambda_x_best(i,:);                % once found, assign the parameters of the agent that is better approximating the i-th subproblem to the parameters of the agent which will be tackling this active subproblem
                f(id_pop_act_subpr(i),:)=lambda_f_best(i,:);                % and do the same for the function values
                id_local=id_local(id_local~=id_pop_act_subpr(i));           % REMOVE this agent from the list of agents performing LOCAL actions (will change when also the IDs of the agents computing SOCIAL actions will change, at subproblems update)
                                
            end
            
        else                                                                % if archsize has only 1 item
            
            id_pop_act_subpr=randperm(n_social);                            % shuffle the vector which contains the IDs of the agents which will perform LOCAL actions
            id_local=setdiff(1:options.popsize,id_pop_act_subpr);           % and REMOVE those IDs from the list of IDs of the actors which will perform SOCIAL actions
            
        end
        
    else                                                                    % if the number of agents performing LOCAL actions is less than the number of objectives (we cannot span the whole space)
        
        id_pop_act_subpr=[];                                                % the vector containing the IDs of the agents performing LOCAL actions on the subproblems will be VOID
        act_subpr=[];                                                       % no subproblem associated
        id_local=1:options.popsize;                                        % ALL the agents will perform SOCIAL actions
        n_social=0;                                                          % and no LOCAL agents will exist anymore
        
    end
    
end

%%  EXPLORE PARAMETERS SETTING

explore_params.F = options.F;
explore_params.CR = options.CR;
explore_params.vlb = vlb;
explore_params.vub = vub;
explore_params.coord_ratio = options.coord_ratio;
explore_params.contr_ratio = options.contr_ratio;
explore_params.rhoini = options.rhoini;
explore_params.cp = options.cp;
explore_params.cpat = options.cpat;
explore_params.ncon = ncon;
explore_params.DE_strategy = options.explore_DE_strategy;
explore_params.func = func;
explore_params.arg = varargin;
explore_params.v = options.v;

%%  SOCIAL PARAMETERS SETTING

social_params.F = options.F;
social_params.CR = options.CR;
social_params.vlb = vlb;
social_params.vub = vub;
social_params.func = func;
social_params.arg = varargin;
social_params.cp = options.cp;
social_params.T = n_social;
social_params.DE_strategy = options.social_DE_strategy;

%%  DRAWING INITIALIZATION

if options.draw_flag>=1                     %if draw_flag
    %tmp_str=func2str(func);        %unused...
    %if draw_flag>1                  %is at least 2, will also save a movie
    %    aviobj = avifile('example.avi','compression','None','quality',25); % needs to be changed to the new VideoWriter interface
    %end
    fig=figure;                     %if is only 1, only plot
    %tmp_igd_ax=1;                  %unused...
end



%%  MAIN LOOP

iter=0;

while (nfeval<options.maxnfeval)
        
    iter=iter+1;
    
    if options.explore_all
        expl_agents=1:options.popsize;                                           % expl_agents specifies the id of the agents which will perform local search. If all pop, all of them
    else
        expl_agents=id_local;                                                   % if not, only "social ones" (this is a contraddiction in terms!!!)
    end
    
    %% INDIVIDUALISTIC MOVES, SELECTION AND ARCHIVING
    
    [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho]=explore(expl_agents,memory,x,v,f,cid,nfeval,rho,explore_params);
        
    z=min([z;ftrial;discarded.f],[],1);                                     % update the min of each function (i.e, the columnwise min).
    
    for i=1:length(expl_agents)                                             % for all agents performing local search
        
        if (maxC(i)<=0)                                                     % if current agent is in feasible region

            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
        
                v(i,:)=vtrial(i,:);                                         % store it's velocity
                x(i,:)=xtrial(i,:);                                         % position
                f(i,:)=ftrial(i,:);                                         % objective function value
                cid(i)=maxC(i);                                             % and max constraint violation
                
            end
            
            if i<=size(discarded.f,2)                                       % if discarded is longer than x_trial, f_trial and maxC, so that we can compare it with them
                
                f_tmp=[ftrial(i,:);discarded.f(i,:)];                       % create a temp f vector, pre-appending ftrial to discarded(i)
                x_tmp=[xtrial(i,:);discarded.x(i,:)];                       % idem with position vector
                c_tmp=[maxC(i);discarded.c(i)];                             % constraint violation
                
            else
                
                f_tmp=ftrial(i,:);                                          % create a temp f vector, pre-appending ftrial to discarded(i)
                x_tmp=xtrial(i,:);                                          % idem with position vector
                c_tmp=maxC(i);                                              % constraint violation
                
            end
            
            domtmp=dominance(f_tmp,0);                                      % compute dominance on f_tmp
            x_tmp=x_tmp(domtmp==0,:);                                       % and resize x_tmp keeping only non-dominated elements of domtmp
            f_tmp=f_tmp(domtmp==0,:);                                       % idem with f
            c_tmp=c_tmp(domtmp==0,:);                                       % and c
            
            [domA,domB]=dominant(memory(:,lx+1:lx+mfit),f_tmp);             % compute dominance of memory wrt f_tmp, and viceversa
            memory=memory(domA==0,:);                                       % and keep in archive only non-dominated ones
            
            if sum(domB==0)>0                                               % if there are non-dominated elements in f wrt to memory
                
                if sum(domA==0)<tmp_max_arch                                % and there is still space in archive
                    
                    memory=[memory;x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)]; % append non dominated elements to archive
                    
                else                                                        % if there's no more space in archive
                    
                    memory=spars(memory,x_tmp(domB==0,:),f_tmp(domB==0,:), c_tmp(domB==0), z); % substitute elements in archive with these ones if doing so would reduce archive clustering in criteria space
                    
                end
                
            end
            
        elseif (maxC(i)>0)&&(maxC(i)<cid(i))                                % if current agent is not in feasible region but it's position improves previous constraint violation
            
            v(i,:)=vtrial(i,:);                                             % store it's velocity
            x(i,:)=xtrial(i,:);                                             % position
            f(i,:)=ftrial(i,:);                                             % objective function values
            cid(i)=maxC(i);                                                 % and max constraint violation
            
            if i<=size(discarded.f,2)
                
                f_tmp=[ftrial(i,:);discarded.f(i,:)];                           % create a temp f vector, pre-appending ftrial to discarded(i)
                x_tmp=[xtrial(i,:);discarded.x(i,:)];                           % idem with position vector
                c_tmp=[maxC(i);discarded.c(i)];                                 % constraint violation
                
            else
                
                f_tmp=ftrial(i,:);                           % create a temp f vector, pre-appending ftrial to discarded(i)
                x_tmp=xtrial(i,:);                           % idem with position vector
                c_tmp=maxC(i);                                 % constraint violation
                
            end
           
            domtmp=dominance(c_tmp,0);                                      % compute the dominance IN CONSTRAINT VIOLATION SENSE of tmp agents
            x_tmp=x_tmp(domtmp==0,:);                                       % and keep only non-dominated elements
            f_tmp=f_tmp(domtmp==0,:);                                       %
            c_tmp=c_tmp(domtmp==0,:);                                       %
            
            [domA,domB]=dominant(memory(:,lx+mfit+2),c_tmp);                % finally, compute dominance IN CONSTRAINT VIOLATION SENSE of tmp agents wrt agents in archive, and viceversa
            memory=memory(domA==0,:);                                       % keep only non dominated elements of memory in memory
            
            if sum(domB==0)>0                                               % if there are non-dominated elements IN CONSTRAINT VIOLATION SENSE in tmp vector
            
                if sum(domA==0)<tmp_max_arch                                % if there's enough space in archive
                
                    memory=[memory;x_tmp(domB==0,:) f_tmp(domB==0,:) zeros(sum(domB==0),1) c_tmp(domB==0)]; % simply append this agent in memory
                    
                else                                                        % if there's no more space in archive
                    
                    memory=spars(memory,x_tmp(domB==0,:),f_tmp(domB==0,:), c_tmp(domB==0), z); % substitute current agents to agents in memory if doing so would reduce clustering in memory
                    
                end
                
            end
            
        end
        
    end
    
    %% SOCIAL MOVES, SELECTION AND ARCHIVING
    
    %delta=max(memory(:,lx+1:lx+mfit),[],1)-min(memory(:,lx+1:lx+mfit),[],1);    % excursion of objective function values
    
    for i=1:n_social                                                         % for each agent tackling a sub-problem
        
        [xsamp,fsamp,maxCsamp,nfeval]=social(x(id_pop_act_subpr(i),:),f(id_pop_act_subpr(i),:),memory,x,cid(id_pop_act_subpr(i)),nfeval,social_params); % perform social actions for i-th subproblem
        
        
        z=min([z;fsamp],[],1);                                              % update minimas
        
        if (maxCsamp<=0)                                                    % if this social action is in feasible region
            
            if ((g_fun(fsamp,lambda(act_subpr(i),:),z)<g_fun(f(id_pop_act_subpr(i),:),lambda(act_subpr(i),:),z))) % if sampled position for agent solving i-th subproblem is good (i.e it's gfun is better than the gfun for the old agent)
                
                if options.v
                
                    v(id_pop_act_subpr(i),:)=xsamp-x(id_pop_act_subpr(i),:);
                    
                end
                
                f(id_pop_act_subpr(i),:)=fsamp;                             % f(agent(active_subproblem)) becomes f of this agent
                x(id_pop_act_subpr(i),:)=xsamp;                             % idem for position
                cid(id_pop_act_subpr(i))=maxCsamp;                          % and constraint violation
                
            end
            
            % Adds sample to global archive
            [domA,domB]=dominant(memory(:,lx+1:lx+mfit),fsamp);             % compute dominance OF OBJECTIVE FUNCTION VALUES of agents in memory wrt new agent and viceversa
            memory=memory(domA==0,:);                                       % restrict memory to non-dominated elements only
            
            if any(domB==0)                                                 % if there's at least one non-dominated element in B (i.e if sampl agent is non-dominated by memory)
                
                if sum(domA==0)<tmp_max_arch                                % and if number of elements in archive is less than allowed (temporarily increased) max
                    
                    memory=[memory;xsamp fsamp domB maxCsamp];              % add this agent to the archive
                    archsize=sum(domA==0)+1;                                % and record increased archive size
                    
                else                                                        % if there's no more space in archive
                    
                    memory=spars(memory,xsamp,fsamp,maxCsamp, z);           % substitute current agent with one in archive if doing so would reduce clustering of agents in archive
                    
                end
                
            end
            
        elseif (maxCsamp>0)&&(maxCsamp<cid(id_pop_act_subpr(i)))            % if this social action is not in feasible region but it improves previous constraint violation
            
            if options.v           
                
                v(id_pop_act_subpr(i),:)=xsamp-x(id_pop_act_subpr(i),:);
                
            end
            
            f(id_pop_act_subpr(i),:)=fsamp;                                 % fsamp becomes f(active_problem) 
            x(id_pop_act_subpr(i),:)=xsamp;                                 % idem with it's position
            cid(id_pop_act_subpr(i))=maxCsamp;                              % and max costraint violation
            
            % Adds sample to global archive
            [domA,domB]=dominant(memory(:,lx+mfit+2),maxCsamp);             % compute dominance OF CONSTRAINT VIOLATION of agents in memory wrt new agent and viceversa 
            memory=memory(domA==0,:);                                       % restrict memory to non dominated agents only (in CONSTRAINT VIOLATION SENSE) 
            
            if any(domB==0)                                                 % if there's at least one non-dominated element in B (i.e if sampl agent is non-dominated by memory)
                
                if sum(domA==0)<tmp_max_arch                                % if there's space in memory
                    
                    memory=[memory;xsamp fsamp domB maxCsamp];              % add current agent to memory
                    archsize=sum(domA==0)+1;                                % and increase current number of agents in memory
                    
                else                                                        % if there's no more space in memory
                    
                    memory=spars(memory,xsamp,fsamp,maxCsamp, z);           % substitute current agent with one in archive if doing so would reduce clustering of agents in archive
                    
                end
                
            end
            
        end
        
    end
    
    sel=memory(:,lx+mfit+1)==0;                                             % select agents with dominance index 0
    memory=memory(sel,:);                                                   % restrict archive to those elements only
    archsize=sum(sel);                                                      % and update archive size
    
    %% ARCHIVE SHRINKING
        
    if (archsize>tmp_max_arch)&&first_shrink                                % if archive is bigger than max allowed within loop, and it's requested to performi a first time shrinking
        
        memory=arch_shrk2(memory,lx,mfit,tmp_max_arch);                      % srhink memory
        archsize=tmp_max_arch;                                              % update size
        first_shrink=0;                                                     % and never do this in-loop shrinking again
        
    end
    
    %% SUBPROBLEMS UPDATE
    
    if n_social&&(mod(iter,supbr_upd_freq)==0)                               % if there are active subproblems and they've been active for the last subpr_upd_freq iterations
        
        [act_subpr,id_pop_act_subpr,pigr,lambda_f_best]=upd_act_subpr(pigr,lambda,lambda_f_best,x,f,memory,z,n_social); % update subproblems
        id_local=setdiff(1:options.popsize,id_pop_act_subpr);               % set the IDs of agents which are performing local actions (i.e differentiate them from the agents which will tacle subproblems through social actions)
        
    end
    
    %% OUTPUT ON SCREEN
    
    if options.draw_flag>=1
        
        hold off
        plot3(memory(:,lx+1),memory(:,lx+2),memory(:,lx+3),'k.')
        hold on
        if (size(id_pop_act_subpr,1)>2)
            plot3(f(id_pop_act_subpr(1:3),1),f(id_pop_act_subpr(1:3),2),f(id_pop_act_subpr(1:3),3),'ro')
            if(size(id_pop_act_subpr,1)>3)
                plot3(f(id_pop_act_subpr(4:end),1),f(id_pop_act_subpr(4:end),2),f(id_pop_act_subpr(4:end),3),'r*')
            end
        end
        drawnow
        
        subplot(2,2,2)
        plot(0:1,options.CR/max(options.CR))
        drawnow
        subplot(2,2,3)
        
        [~,id_pigr]=sort(-pigr);
        
        switch mfit
            
            case 2
                
                plot(lambda(:,1),lambda(:,2),'.')
                hold on
                plot(lambda(act_subpr,1),lambda(act_subpr,2),'ro')
                plot(lambda(id_pigr(1:n_social),1),lambda(id_pigr(1:n_social),2),'ko')
                
            case 3
                
                plot3(lambda(:,1),lambda(:,2),lambda(:,3),'.')
                hold on
                plot3(lambda(act_subpr,1),lambda(act_subpr,2),lambda(act_subpr,3),'r*')
                plot3(lambda(id_pigr(1:n_social),1),lambda(id_pigr(1:n_social),2),lambda(id_pigr(1:n_social),3),'ko')
                
        end
        drawnow
        %         figure(1)
        subplot(2,2,4)
        [~,kkk]=sort(lambda(:,1));
        hold off
        switch mfit
            case 2
                %                     hold(AX(1),'off')
                [AX,H1,H2]=plotyy(1:n_lambda,pigr(kkk),1:options.popsize,rho,'plot','semilogy');
                set(H1,'LineStyle','-.','color','b')
                set(H2,'LineStyle','-.','color','g')
                %                     hold(AX(1),'on')
                set(AX(1),'YLim',[0 1])
                %set(AX(2),'YLim',[tolconv rhoini],'YGrid','on')
                %axis(AX(2),[0 popsize tolconv 0.5])
                %                     hold(AX(1),'on')
                %                     drawnow
            case 3
                betat=asin(lambda(:,3));
                alphat=atan(lambda(:,2)./lambda(:,1));
                plot3(alphat,betat,pigr,'.')
                axis([0 pi/2 0 pi/2 0 1])
                hold on
        end
        %          figure(1)
        subplot(2,2,1)
        switch mfit
            case 2
                plot(memory(:,lx+1),memory(:,lx+2),'k.',f(:,1),f(:,2),'ro')
                hold on
                plot(f(id_pop_act_subpr(1:2),1),f(id_pop_act_subpr(1:2),2),'gd')
                plot(f(id_pop_act_subpr(3:end),1),f(id_pop_act_subpr(3:end),2),'r*')
            case 3
                plot3(memory(:,lx+1),memory(:,lx+2),memory(:,lx+3),'k.')
                hold on
                plot3(f(id_pop_act_subpr(1:3),1),f(id_pop_act_subpr(1:3),2),f(id_pop_act_subpr(1:3),3),'ro')
                plot3(f(id_pop_act_subpr(4:end),1),f(id_pop_act_subpr(4:end),2),f(id_pop_act_subpr(4:end),3),'r*')
        end
        drawnow
        hold off
        if options.draw_flag>1
            FF = getframe(fig);
            aviobj = addframe(aviobj,FF);
        end
    end
    
    if (mod(iter,50)==0)&&~isempty(filename)
        save([filename '_tmp'],'memory','x','f','iter')
    end
    
    if options.draw_flag>0
        fprintf('Iter: %d - feval: %d/%d.',iter,nfeval,options.maxnfeval)
        fprintf('\tArch. size: %d.',size(memory,1))
        fprintf('\tNon-dominated pop.: %d/%d.',sum(dominance(f,0)==0),options.popsize)
        fprintf('\n')
    end
    
end

%% ARCHIVE SHRINKING BEFORE OUTPUT

if archsize>options.max_arch
    
    memory=arch_shrk2(memory,lx,mfit,options.max_arch);
    
end

%% PLOTTING, IF REQUIRED

if options.draw_flag>=1
    
    %    figure(1)
    subplot(2,2,1)
    hold off
    switch mfit
        
        case 2
            
            plot(memory(:,lx+1),memory(:,lx+2),'k.')
            hold on
            
        case 3
            
            plot3(memory(:,lx+1),memory(:,lx+2),memory(:,lx+3),'k.')
            hold on
    end
    
    drawnow
    
    if options.draw_flag>1
        
        FF = getframe(fig);
        aviobj = addframe(aviobj,FF);
        aviobj=close(aviobj);
        
    end
end

return
