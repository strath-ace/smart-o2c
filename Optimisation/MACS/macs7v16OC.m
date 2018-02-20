function [memory,nfeval,ener,history]=macs7v16OC(func,memory,vlb,vub,options,filename,fileload,varargin)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Multi Agent Collaborative Search

%  memory=macs7v16OC(func,vlb,vub,options,filename,fileload,varargin)
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
%                   'DE_strategy' Strategy to use for DE (ie pull towards
%                               best element or random one, default 'best')
%
%  OUTPUT
%           memory :  matrix containing the archived solutions and the
%                       associated value of the cost function
%                       memory=[x f]
%           nfeval   :  number of total function evaluations

% CHANGE LOG
% Created by: Federico Zuiani 2011
% Revised, cleaned and optimized by: Lorenzo A. Ricciardi 2015-2018

%%  MACS PARAMETERS DEFAULT VALUES

default.maxnfeval=5000;
default.popsize=10;
default.min_popsize=1;
default.rhoini=1;
default.F=0.9;
default.Fmax=2;
default.CR=1;
default.p_social=1;
default.max_arch=10;
default.coord_ratio=1;
default.contr_ratio=0.5;
default.draw_flag=0;
default.cp=0;
default.MBHflag=0;
default.mbh_options = '';
default.timeout_single_level_opt = inf;
default.explore_DE_strategy='best';
default.social_DE_strategy='DE/current-to-rand/1';
default.v=0;
default.dyn_pat_search = 1;
default.upd_subproblems = 1;
default.max_rho_contr = 5;
default.pat_search_strategy = 'standard';
default.optimal_control = 0;
default.bilevel = 0;
default.vars_to_opt = ones(length(vlb),1);
default.prob_subdiv = 3;
default.cent_prob = 1;

%%  MACS PARAMETERS VALIDATION AND DEFAULTING

vlb = (vlb(:))';        % vlb must be a row vector
vub = (vub(:))';        % vub must be a row vector

if isfield(options,'maxnfeval')
    if ((options.maxnfeval<=0) || mod(options.maxnfeval>0,1)>0) %~isPositiveIntegerValuedNumeric(options.maxnfeval)
        error('maxnfeval must be a positive integer!');
    end
else
    warning(['maxnfeval not supplied, using default value ' num2str(default.maxnfeval)]);
    options.maxnfeval = default.maxnfeval;
end

if isfield(options,'popsize')
    if ((options.popsize<=0) || mod(options.popsize>0,1)>0) %~isPositiveIntegerValuedNumeric(options.popsize)
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
    if ((options.max_arch<=0) || mod(options.max_arch>0,1)>0) %~isPositiveIntegerValuedNumeric(options.max_arch)
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
    if ((options.draw_flag<0) || mod(options.draw_flag>0,1)>0) %~isPositiveIntegerValuedNumeric(options.draw_flag+1) %0 is not counted as postive, this should do the trick
        error('draw_flag must be a non negative integer!');
    end
else
    warning(['draw_flag not supplied, using default value ' num2str(default.draw_flag)]);
    options.draw_flag = default.draw_flag;
end

if isfield(options,'cp')
    if ((options.cp<0) || (options.cp<0) ||  mod(options.cp>0,1)>0) %~isPositiveIntegerValuedNumeric(options.cp+1) %0 is not counted as postive, this should do the trick
        error('cp must be a 0 or 1!');
    end
else
    warning(['cp not supplied, using default value ' num2str(default.cp)]);
    options.cp = default.cp;
end

if isfield(options,'MBHflag')
    if ((options.MBHflag<0) || mod(options.MBHflag>0,1)>0) %~isPositiveIntegerValuedNumeric(options.MBHflag+1) %0 is not counted as postive, this should do the trick
        error('MBHflag must be a non negative integer');
    end
    if options.MBHflag>0
        if ~isfield(options,'mbh_options')
            error('MBHflag>1 but no mbh_options specified');
        end
        if ~isfield(options,'smooth_scal_constr_fun')
            error('MBHflag>1 but no smooth_scal_constr_fun specified');
        end
        if ~isfield(options,'timeout_single_level_opt')
            warning(['MBHflag>1 but no timeout_single_level_opt specified, using default value ' num2str(default.timeout_single_level_opt)]);
            options.timeout_single_level_opt = default.timeout_single_level_opt;
        end
    else
        options.mbh_options = '';
        options.smooth_scal_constr_fun  = '';
        options.timeout_single_level_opt = default.timeout_single_level_opt;
    end
else
    warning(['MBHflag not supplied, using default value ' num2str(default.MBHflag)]);
    options.MBHflag = default.MBHflag;
    options.mbh_options = '';
    options.smooth_scal_constr_fun = '';
    options.timeout_single_level_opt = default.timeout_single_level_opt;
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

if isfield(options,'v')
    if ((options.v<0) || (options.v>1) || mod(options.v>0,1)>0) %~isPositiveIntegerValuedNumeric(options.v+1) %0 is not counted as postive, this should do the trick
        error('v must be a 0 or 1!');
    end
else
    warning(['v not supplied, using default value ' num2str(default.v)])
    options.v = default.v;
end

if isfield(options,'dyn_pat_search')
    if ((options.dyn_pat_search<0) || (options.dyn_pat_search>1) || mod(options.dyn_pat_search>0,1)>0) % ~isPositiveIntegerValuedNumeric(options.dyn_pat_search+1) %0 is not counted as postive, this should do the trick
        error('dyn_pat_searchv must be a 0 or 1!');
    else
        if options.dyn_pat_search>1
            error('dyn_pat_search must be a 0 or 1!');
        end
    end
else
    warning(['dyn_pat_search not supplied, using default value ' num2str(default.dyn_pat_search)])
    options.dyn_pat_search = default.dyn_pat_search;
end

if isfield(options,'upd_subproblems')
    if ((options.upd_subproblems<0) || (options.upd_subproblems>1) || mod(options.upd_subproblems>0,1)>0) % ~isPositiveIntegerValuedNumeric(options.upd_subproblems+1) %0 is not counted as postive, this should do the trick
        error('upd_subproblems must be a 0 or 1!');
    else
        if options.upd_subproblems>1
            error('upd_subproblems must be a 0 or 1!');
        end
    end
else
    warning(['upd_subproblems not supplied, using default value ' num2str(default.upd_subproblems)])
    options.upd_subproblems = default.upd_subproblems;
end

if isfield(options,'max_rho_contr')
    if ((options.max_rho_contr<=0) || mod(options.max_rho_contr>0,1)>0) % ~isPositiveIntegerValuedNumeric(options.max_rho_contr+1)
        error('max_rho_contr must be a non negative integer!');
    end
else
    warning(['max_rho_contr not supplied, using default value ' num2str(default.max_rho_contr)]);
    options.max_rho_contr = default.max_rho_contr;
end

if isfield(options,'pat_search_strategy')
    if ~(strcmp(options.pat_search_strategy,'standard') || strcmp(options.pat_search_strategy,'tracking') || strcmp(options.pat_search_strategy,'random') || strcmp(options.pat_search_strategy,'MADS') || strcmp(options.pat_search_strategy,'none'))
        error('Unrecognized pat_search_strategy: must be either standard, tracking or none')
    end
    if strcmp(options.pat_search_strategy,'MADS')
        %if options.dyn_pat_search
        %    warning('MADS and dynamic pattern search are not compatible. Setting dyn_pat_search to 0');
        %    options.dyn_pat_search = 0;
        %end
    end
else
    warning(['pat_search_strategy not supplied, using defauld value ' default.pat_search_strategy])
    options.pat_search_strategy = default.pat_search_strategy;
end

if isfield(options,'optimal_control')
    if ((options.optimal_control<0) || (options.optimal_control>1) || mod(options.optimal_control>0,1)>0)  % ~isPositiveIntegerValuedNumeric(options.optimal_control+1) %0 is not counted as postive, this should do the trick
        error('optimal_control must be a 0 or 1!');
    else
        if options.optimal_control>1
            error('optimal_control must be a 0 or 1!');
        end
    end
else
    warning(['optimal_control not supplied, using default value ' num2str(default.optimal_control)])
    options.optimal_control = default.optimal_control;
end

if isfield(options,'bilevel')
    if ((options.bilevel<0) || (options.bilevel>1) || mod(options.bilevel>0,1)>0) %~isPositiveIntegerValuedNumeric(options.bilevel+1) %0 is not counted as postive, this should do the trick
        error('bilevel must be a 0 or 1!');
    else
        if options.bilevel>1
            error('bilevel must be a 0 or 1!');
        end
    end
else
    warning(['bilevel not supplied, using default value ' num2str(default.bilevel)])
    options.bilevel = default.bilevel;
end

if isfield(options,'vars_to_opt')
    if (max(size(options.vars_to_opt))~=max(size(vub))) || (min(size(options.vars_to_opt))~=min(size(vub)))
        error('vars_to_opt must have the same number of elements as vlb and vub') ;
    end
    if (any(options.vars_to_opt<0) || any(options.vars_to_opt>1) || any(mod(options.vars_to_opt>0,1)>0)) % ~isPositiveIntegerValuedNumeric(options.vars_to_opt+1) %0 is not counted as postive, this should do the trick
        error('all vars_to_opt entries must be a 0 or 1!');
    end
    if all(options.vars_to_opt==0)
        error('all vars_to_opt are zero: no oprimisation possible');
    end
else
    warning('vars_to_opt not supplied, using defauld value (all)')
    options.vars_to_opt = default.vars_to_opt;
end

options.vars_to_opt = logical(options.vars_to_opt);

refine = 0;

if options.MBHflag>0
    
    if ~isfield(options,'smooth_scal_constr_fun')
        
        error('No smooth_scal_constr_fun field found');
        
    end
    
    if ~isfield(options,'refine_freq')
        
        error('MBHflag specified but no refine_freq field found');
        
    else
        
        if ((options.refine_freq<=0) || mod(options.refine_freq>0,1)>0) %~isPositiveIntegerValuedNumeric(options.refine_freq)
            
            error('refine_freq be a positive integer!');
            
        else
            
            refine = 1;
            
        end
        
    end
    
end

if options.optimal_control==1
    
    if isfield(options,'oc')
        
        if ~isfield(options.oc,'problem')
            
            error('oc.problem must be supplied');
            
        end
        
    else
        
        error('optimal_control flag is 1, but no oc structure is supplied');
        
    end
    
else
    
    scales = max([abs(vub);abs(vlb);abs(vub-vlb)],[],1);
    
end

if isfield(options,'prob_subdiv')
    if ((options.prob_subdiv<=0) || mod(options.prob_subdiv>0,1)>0) %~isPositiveIntegerValuedNumeric(options.prob_subdiv+1) %0 is not counted as postive, this should do the trick
        error('prob_subdiv must be a non negative integer');
    end
else
    warning(['prob_subdiv not supplied, using default value ' num2str(default.prob_subdiv)]);
    options.prob_subdiv = default.prob_subdiv;
end

if isfield(options,'cent_prob')
    if ((options.cent_prob<0) || (options.cent_prob>1) || mod(options.cent_prob>0,1)>0) %~isPositiveIntegerValuedNumeric(options.cent_prob+1) %0 is not counted as postive, this should do the trick
        error('cent_prob must be 0 or 1');
    end
else
    warning(['cent_prob not supplied, using default value ' num2str(default.cent_prob)]);
    options.cent_prob = default.cent_prob;
end

%%  MACS SUBPARAMETERS AUTOSETTING

n_social        = round(options.p_social*options.popsize);                  % number of agents performing social actions
index = 1:length(vlb);
options.id_vars_to_opt = index(options.vars_to_opt);

if options.upd_subproblems
    
    supbr_upd_freq = n_social;                                                  % number of iterations after which new subproblems are chosen
    
else
    
    supbr_upd_freq = inf;                                                  % number of iterations after which new subproblems are chosen
    
end

lx = length(vlb);                                               % number of parameters of objective functions
dd = [];                                                                    % matrix of pairwise energies, needed for storage
energy = 0;                                                                 % archive energy, needed for storage
ener = [];
ener2 = [];

%%  MEMORY INITIALIZATION

nfeval=0;

if isempty(memory)
    
    %% CHOOSE RANDOM PARAMETERS AND SET UP PROBLEM DIMENSIONALITY AND CONSTRAINTS
    
    if options.optimal_control
        
        % Random initial guess for all
        
        % Generate random controls and static variables
        
        x=lhsu(vlb,vub,options.popsize,options.id_vars_to_opt);                                        % latin hypercube sampling over the whole domain
        
        % Generate a near feasible initial condition for the given static
        % and control variables
        
        for i = 1:options.popsize
            
            fprintf('Generating feasible initial guess for agent %d...\n',i);
            x(i,:) = generate_guess(options.oc.problem,0,x(i,:));
            
        end
        
    else
        
        x=lhsu(vlb,vub,options.popsize,options.id_vars_to_opt);                                        % latin hypercube sampling over the whole domain
        
    end
    
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
    
    c = zeros(options.popsize,ncon);                                        % initialize the matrix of constraint violations
    mins = realmax*ones(1,mfit);
    maxs = zeros(1,mfit);
    
    %% CREATION OF INITIAL POPULATION
        
    f=zeros(options.popsize,mfit);                                      % initialize a matrix of zeros as f
    
    for i=1:options.popsize                                             % and populate it:
        
        fprintf(['Initialising agent ',num2str(i),'\n']);
        
        y=x(i,:);                                                       % create a temporary y
        
        if options.cp==0                                                % if the problem is unconstrained
            
            if options.bilevel == 1
                
                [f(i,:),xcorr]=func(y,varargin{:});                         % evaluate f
                x(i,:) = xcorr;                                             % substitute initial random x with FEASIBLE x
                
            else
                
                f(i,:)=func(y,varargin{:});                         % evaluate f
                
            end
            
            maxC(i)=0;                                                  % and say all the constraints are satisfied
            
        else                                                            % if the problem is constrained
            
            [f(i,:),c(i,:)]=func(y,varargin{:});                        % evaluate f
            maxC(i)=max(c(i,:));                                        % and say maximum violation is max of c
            
        end
        
        nfeval=nfeval+1;                                                % then update nfeval
        
    end
    
    dom=dominance(f,0);                                                 % once you have f, evaluate dominances. BEWARE, THIS SHOULD ALSO BE DONE FOR LOADED PREVIOUS EVALUATIONS!!!
    cid=maxC;
    qq = [x(dom==0,:) f(dom==0,:) dom(dom==0) cid(dom==0,:)];
    
     if size(qq,1)<=options.max_arch
    
        [memory,dd,energy,ener2,mins,maxs]=arch_shrk6([],dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        
    else
        
        [mm,dd,energy,ener2,mins,maxs]=arch_shrk6([],dd,qq(1:options.max_arch,:),energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(mm,dd,qq(options.max_arch+1:end,:),energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        
    end
      
else
    
    %% LOAD MEMORY AND CHECK CONSISTENCY WITH THE VARIABLES
    
    [n_rows,n_cols] = size(memory);
    
    n_obj_trial = n_cols-lx-2;
    
    if n_obj_trial<1        % memory should at least have vlb+3 columns (if the objective function retuns a scalar)
        
        error('The archive memory should at least have 3 columns more than the length of the vector of lower bounds for the variables')
        
    else
        
        x = memory(:,1:lx);
        f = memory(:,lx+1:lx+n_obj_trial);
        cid = memory(:,end);
        
    end
    
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
    
    if mfit ~= n_obj_trial
        
        error('The number of objectives deducted from the archive is not the same as the one deducted from calling the objective function')
        
    end
    
    mins = realmax*ones(1,mfit);
    maxs = zeros(1,mfit);
    
    %% IF SUPPLIED ARCHIVE IS SMALLER THAN WORKING POPULATION, RANDOMLY CREATE NEW AGENTS
    
    if n_rows<options.popsize
        
        warning('Supplied population is smaller than the requested number of agents, creating the remaining ones through Latin Hypercube Sampling');
        
        fadd=zeros(options.popsize-n_rows,mfit);                                      % initialize a matrix of zeros as f
        cadd = zeros(options.popsize-n_rows,ncon);
        maxCadd = zeros(options.popsize-n_rows,1);
        
        if options.optimal_control
            
            % Random initial guess for all
            
            % Generate random controls and static variables
            
            xadd=lhsu(vlb,vub,options.popsize-n_rows,options.id_vars_to_opt);                                        % latin hypercube sampling over the whole domain
            
            % Generate a near feasible initial condition for the given static
            % and control variables
            
            for i = 1:options.popsize-n_rows
                
                fprintf('Generating feasible initial guess for agent %d...\n',nrows+i);
                xadd(i,:) = generate_guess(options.oc.problem,0,xadd(i,:));
                
            end
            
        else
            
            xadd=lhsu(vlb,vub,options.popsize-n_rows,options.id_vars_to_opt);                                        % latin hypercube sampling over the whole domain
            
        end
        
        for i=1:options.popsize-n_rows                                                       % and populate it:
            
            fprintf(['Initialising agent ',num2str(n_rows+i),'\n']);
            
            y=xadd(i,:);                                                       % create a temporary y
            
            if options.cp==0                                                % if the problem is unconstrained
                
                if options.bilevel == 1
                    
                    [fadd(i,:),xcorr]=func(y,varargin{:});                         % evaluate f
                    xadd(i,:) = xcorr;                                             % substitute initial random x with FEASIBLE x
                    
                else
                    
                    fadd(i,:)=func(y,varargin{:});                         % evaluate f
                    
                end
                
                maxCadd(i)=0;                                                  % and say all the constraints are satisfied
                
            else                                                            % if the problem is constrained
                
                [fadd(i,:),cadd(i,:)]=func(y,varargin{:});                        % evaluate f
                maxCadd(i)=max(cadd(i,:));                                        % and say maximum violation is max of c
                
            end
            
            nfeval=nfeval+1;                                                % then update nfeval
            
        end
        
        f = [f;fadd];
        x = [x;xadd];
        cid = [cid;maxCadd];
        
    end
    
    %% FILTER DOMINATED SOLUTIONS 
    
    dom=dominance(f,0);                                                 % get nondominated solutions
    qq = [x(dom==0,:) f(dom==0,:) dom(dom==0) cid(dom==0,:)];
    
    if size(qq,1)<=options.max_arch
    
        [memory,dd,energy,ener2,mins,maxs]=arch_shrk6([],dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        
    else
        
        [mm,dd,energy,ener2,mins,maxs]=arch_shrk6([],dd,qq(1:options.max_arch,:),energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(mm,dd,qq(options.max_arch+1:end,:),energy,ener2,mins,maxs,lx,mfit,options.max_arch);
        
    end

    % f and x can now only be >= than popsize. If thery're larger, we need
    % some strategy to cull the number
    
    if size(f,1)>options.popsize
        
        if size(memory,1)>=options.popsize      % if archive has more points than popsize, we get just what we need
    
            [mm,ddd,energy,ener2,mins,maxs]=arch_shrk6([],[],memory(1:options.popsize,:),0,[],mins,maxs,lx,mfit,options.popsize);
            [mm,ddd,energy,ener2,mins,maxs]=arch_shrk6(mm,ddd,memory(options.popsize+1:end,:),energy,ener2,mins,maxs,lx,mfit,options.popsize);

            x = mm(:,1:lx);
            f = mm(:,lx+1:lx+mfit);
            cid = mm(:,end);
            
        else
            
            % we get everything we can from the archive, and select 
            % randomly the remaining ones
            
            xx = memory(:,1:lx);
            ff = memory(:,lx+1:lx+mfit);
            cc = memory(:,end);
            
            fdiff = setdiff(f,ff,'rows');
            
            idchoose = randperm(size(fdiff,1));
            idchoose = idchoose(1:options.popsize-size(ff,1));
            
            f = [ff;f(idchoose,:)];
            x = [xx;x(idchoose,:)];
            cid = [cc,cid(idchoose)];
            
        end
    

    end
    
    
end

delta=max(memory(:,lx+1:lx+mfit),[],1)-min(memory(:,lx+1:lx+mfit),[],1);    % compute delta, the excursion in criteria space for each element of the population(can also be 0!)
z=min(memory(:,lx+1:lx+mfit),[],1);                                   % and z, the min in criteria space NOTE, THIS IS NOT THE TRUE MIN, ONLY THE BEST COMPUTED SO FAR!!!!
zstar=max([memory(:,lx+1:lx+mfit);f],[],1);
archsize=length(memory(:,1));

%% INITIALIZATION OF LAMBDA VECTORS DEFINING SUBPROBLEMS

switch mfit
    
    case 1                                                              % Single Objective
        
        n_lambda=n_social;                                                     % number of subproblems
        lambda=ones(n_social,1);%[1;1];                                                   % at least 2 subproblems are needed for the whole strategy to work
        
    case 2                                                              % Bi-objective
        
        if options.upd_subproblems                                                  % number of iterations after which new subproblems are chosen
            
            n_lambda=mfit*100;                                              % number of subproblems is 100 times the number of objectives, so always 200... that 100 could be a parameter...
            
        else
            
            n_lambda=n_social;
            
        end
        
        alfas=linspace(0,pi/2,n_lambda)';                               % subproblems of a Bi-objective optimization are distributed regularly around a 1/4 circle, so here we compute the angles
        alfas=alfas(2:end-1);                                           % orthogonal subproblems (alpha=0 and pi/2) are excluded from this list, since they will be introduced manually later
        alfas=alfas(randperm(n_lambda-2));                              % shuffle the alphas, to shuffle the subproblems (why is this needed???)
        lambda=[eye(mfit); sin(alfas) cos(alfas)];                      % finally, in lambdas are m orthogonal subproblems, plus all the other ones around a 1/4 circle
        
    otherwise                                                           % Many objectives
        
        lambda = generate_uniform_directions(mfit,options.prob_subdiv,options.cent_prob);
        n_lambda = size(lambda,1);
        
end

%% INITIALIZATION OF PIGR, V AND RHO

pigr=ones(1,n_lambda);              % utility function
v=zeros(options.popsize,lx);                % PROBABLY THIS IS THE VECTOR OF IMPROVEMENT SINCE THE LAST ITERATION
rho(:,1)=options.rhoini*ones(options.popsize,1);    %the radius of the local search
rho(:,2)=zeros(options.popsize,1);          %counter of iterations, after some will reset radius to initial value

%% CHECK WHICH IS THE BEST APPROXIMATION OF THE I-TH SUBPROBLEM

lambda_f_best=zeros(n_lambda,mfit);                                     % initialize lambda_f (same as criteria space dimensionality)
lambda_x_best=zeros(n_lambda,lx);                                       % initialize lambda_x (same as parameter space dimensionality)

for i=1:n_lambda                                                        % for each subproblem
    
    g_tmp = hp_g_fun(memory(:,lx+1:lx+mfit),lambda(i,:),z,zstar);       % evaluate g_fun on all agents for current i-th subproblem
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
            
            %x(id_pop_act_subpr(i),:)=lambda_x_best(i,:);                % once found, assign the parameters of the agent that is better approximating the i-th subproblem to the parameters of the agent which will be tackling this active subproblem
            %f(id_pop_act_subpr(i),:)=lambda_f_best(i,:);                % and do the same for the function values
            id_local=id_local(id_local~=id_pop_act_subpr(i));           % REMOVE this agent from the list of agents performing LOCAL actions (will change when also the IDs of the agents computing SOCIAL actions will change, at subproblems update)
            
        end
        
    else                                                                % if archsize has only 1 item
        
        id_pop_act_subpr=randperm(n_social);                            % shuffle the vector which contains the IDs of the agents which will perform LOCAL actions
        id_local=setdiff(1:options.popsize,id_pop_act_subpr);           % and REMOVE those IDs from the list of IDs of the actors which will perform SOCIAL actions
        
    end
    
    if any(id_pop_act_subpr==0)
        error('Problems with association subproblem-agent');
    end
    
else                                                                    % if the number of agents performing LOCAL actions is less than the number of objectives (we cannot span the whole space)
    
    id_pop_act_subpr=[];                                                % the vector containing the IDs of the agents performing LOCAL actions on the subproblems will be VOID
    act_subpr=[];                                                       % no subproblem associated
    id_local=1:options.popsize;                                        % ALL the agents will perform SOCIAL actions
    n_social=0;                                                          % and no LOCAL agents will exist anymore
    
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
explore_params.DE_strategy = options.explore_DE_strategy;
explore_params.func = func;
explore_params.arg = varargin;
explore_params.v = options.v;
explore_params.max_rho_contr = options.max_rho_contr;
explore_params.MBHflag = options.MBHflag;
explore_params.mbh_options = options.mbh_options;
explore_params.pat_search_strategy= options.pat_search_strategy;
explore_params.optimal_control = options.optimal_control;
explore_params.vars_to_opt = options.vars_to_opt;
explore_params.id_vars_to_opt = options.id_vars_to_opt;
explore_params.id_vars_to_leave = setdiff(1:length(vlb),options.vars_to_opt);
explore_params.foptionsNLP = options.mbh_options;
explore_params.smooth_scal_constr_fun = options.smooth_scal_constr_fun;
explore_params.bilevel = options.bilevel;
explore_params.timeout_single_level_opt = options.timeout_single_level_opt;

if explore_params.optimal_control == 1
    
    explore_params.oc = options.oc;
    
else
    
    explore_params.scales = scales;
    explore_params.norm_vub = vub./scales;
    explore_params.norm_vlb = vlb./scales;
    explore_params.bilevel = options.bilevel;
    explore_params.mbh_options = options.mbh_options;
    explore_params.mbh_func = options.func;
    explore_params.mbh_cfunc = options.cfunc;
    
end

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
social_params.vars_to_opt = options.vars_to_opt;
social_params.id_vars_to_opt = options.id_vars_to_opt;
social_params.optimal_control = options.optimal_control;
social_params.bilevel = options.bilevel;

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

history = [0*ones(size(memory,1),1) memory];

fprintf(['Initialisation finished, ', num2str(size(memory,1)), ' non dominated solutions found... \n']);
fprintf('Starting optimisation loop...\n');

MBH_positions = nan(1,lx+mfit); %positions and lambdas
MADS_dirs = [];
int_loc_opt = 0;

iter=0;
oldmin = min(memory(:,lx+1:lx+mfit),[],1);

for i =1:options.popsize
    
    patdirs(i).avail = 1:lx;
    
end

oldmem = memory;

while nfeval<options.maxnfeval
    
    fold = f;
    
    nfeval_old = nfeval;
    
    tic
    
    newmins = (min(memory(:,lx+1:lx+mfit))<oldmin);
    
    
    % Pattern search coord_ratio dynamically adapted, depending on filling
    % ratio of archive
    % Currently begins scanning along all coordinates and when the archive
    % is full looks for only one coordinate.
    % Suggestion: MIGHT be useful to memorise the coordinates already
    % scanned and always pick for a different one if the agent hasn't moved
    
    if options.dyn_pat_search
        
        fillratio = size(memory,1)/options.max_arch;
        patsearmin = 1/lx;
        explore_params.coord_ratio = 1-(1-patsearmin)*fillratio;
        
    end
    
    iter=iter+1;
    
    %% INDIVIDUALISTIC MOVES, SELECTION AND ARCHIVING
    
    local_only = 0;
    
    explore_params.check_grad = 0;
    
    [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho,patdirs,MBH_positions,MADS_dirs]=explore2(memory,x,v,f,cid,nfeval,lambda,act_subpr,id_pop_act_subpr,z,zstar,rho,patdirs,pigr,MBH_positions,MADS_dirs,local_only,int_loc_opt,explore_params);
    
    oldz = z;
    z=min([z;ftrial;discarded.f],[],1);            % update the min of each function (i.e, the columnwise min).
    zstar=max([memory(:,lx+1:lx+mfit); f],[],1);
    
    if any(z<oldz)
        
        fprintf('INDIVIDUAL MOVES IMPROVED MINIMA\n');
        
        for i=1:mfit
            
            if z(i)<oldz(i)
                
                fprintf('Min of objective %d has changed! Old value = %g, new value = %g \n',i,oldz(i),z(i));
                
            end
            
        end
        
    end
    
    for i=1:options.popsize                                             % for all agents performing local search
        
        if (maxC(i)<=0)                                                     % if current agent is in feasible region
            
            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                
                v(i,:) = vtrial(i,:);                                       % given velocity
                
                % nornalisation of v, to avoid friction
                %                 if norm(v(i,:))>0
                %
                %                     v(i,:) = v(i,:)/norm(v(i,:));
                %
                %                 end
                
                x(i,:)=xtrial(i,:);                                         % position
                f(i,:)=ftrial(i,:);                                         % objective function value
                cid(i)=maxC(i);                                             % and max constraint violation
                patdirs(i).avail = 1:lx;
                
            end
            
        elseif (maxC(i)>0)&&(maxC(i)<cid(i))                                % if current agent is not in feasible region but it's position improves previous constraint violation
            
            v(i,:)=xtrial(i,:)-x(i);                                             % store it's velocity
            x(i,:)=xtrial(i,:);                                             % position
            f(i,:)=ftrial(i,:);                                             % objective function values
            cid(i)=maxC(i);                                                 % and max constraint violation
            
        end
        
    end
    
    % splitting feasible and non feasible solutions
    ffeas = discarded.f(discarded.c<=0,:);
    xfeas = discarded.x(discarded.c<=0,:);
    cfeas = discarded.c(discarded.c<=0,:);
    
    fnfeas = discarded.f(discarded.c>0,:);
    xnfeas = discarded.x(discarded.c>0,:);
    cnfeas = discarded.c(discarded.c>0,:);
    
    if ~isempty(ffeas)
        
        domtmp=dominance(ffeas,0);                                          % compute dominance on ffeas
        x_tmp=xfeas(domtmp==0,:);                                           % and resize x_tmp keeping only non-dominated elements of domtmp
        f_tmp=ffeas(domtmp==0,:);                                           % idem with f
        c_tmp=cfeas(domtmp==0,:);                                           % and c
        
        [domA,domB]=dominant(memory(memory(:,end)<=0,lx+1:lx+mfit),f_tmp);  % compute dominance of memory wrt f_tmp, and viceversa
        
        [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
        
        if sum(domB==0)>0                                                   % if there are non-dominated elements in f wrt to memory
            
            qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
            
            [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
            
        end
        
    end
    
    
    if ~isempty(fnfeas)
        
        domtmp=dominance(cnfeas,0);                                         % compute dominance on ftrial
        x_tmp=xnfeas(domtmp==0,:);                                          % and resize x_tmp keeping only non-dominated elements of domtmp
        f_tmp=fnfeas(domtmp==0,:);                                          % idem with f
        c_tmp=cnfeas(domtmp==0,:);                                          % and c
        
        [domA,domB]=dominant(memory(:,end),c_tmp);                               % compute dominance of memory wrt cf_tmp, and viceversa
        [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
        
        if sum(domB==0)>0                                                       % if there are non-dominated elements in f wrt to memory
            
            qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
            
            [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
            
        end
        
    end
    
    %% SOCIAL MOVES, SELECTION AND ARCHIVING
    
    xsoc = zeros(n_social,lx);
    fsoc = zeros(n_social,mfit);
    csoc = zeros(n_social,ncon);
    
    fprintf(['Performing social actions: ']);
    
    for i=1:n_social                                                        % for each agent tackling a sub-problem
        
        fprintf([num2str(i),', ']);
        
        [xsamp,fsamp,maxCsamp,nfeval]=social2(x(id_pop_act_subpr(i),:),f(id_pop_act_subpr(i),:),memory,x,cid(id_pop_act_subpr(i)),nfeval,energy,ener2,social_params); % perform social actions for i-th subproblem
        
        xsoc(i,:) = xsamp;
        fsoc(i,:) = fsamp;
        csoc(i,:) = maxCsamp;
        
        % Old style social move. If a sample generated through social DE is
        % "good" (i.e. minimises the current agent's subproblem) then it
        % moves in this position. Not the best choice because a better
        % position might be found by a social or local action of another
        % agent, new implementation is much better.
        
        %                         if (maxCsamp<=0)                                                    % if this social action is in feasible region
        %
        %                             if ((g_fun(fsamp,lambda(act_subpr(i),:),z)<g_fun(f(id_pop_act_subpr(i),:),lambda(act_subpr(i),:),z))) % if sampled position for agent solving i-th subproblem is good (i.e it's gfun is better than the gfun for the old agent)
        %
        %                                 if options.v
        %
        %                                     v(id_pop_act_subpr(i),:)=xsamp-x(id_pop_act_subpr(i),:);
        %
        %                                 end
        %
        %                                 f(id_pop_act_subpr(i),:)=fsamp;                             % f(agent(active_subproblem)) becomes f of this agent
        %                                 x(id_pop_act_subpr(i),:)=xsamp;                             % idem for position
        %                                 cid(id_pop_act_subpr(i))=maxCsamp;                          % and constraint violation
        %
        %                             end
        %
        %                         elseif (maxCsamp>0)&&(maxCsamp<cid(id_pop_act_subpr(i)))            % if this social action is not in feasible region but it improves previous constraint violation
        %
        %                             if options.v
        %
        %                                 v(id_pop_act_subpr(i),:)=xsamp-x(id_pop_act_subpr(i),:);
        %
        %                             end
        %
        %                             f(id_pop_act_subpr(i),:)=fsamp;                                 % fsamp becomes f(active_problem)
        %                             x(id_pop_act_subpr(i),:)=xsamp;                                 % idem with it's position
        %                             cid(id_pop_act_subpr(i))=maxCsamp;                              % and max costraint violation
        %
        %                         end
        
    end
    
    fprintf('\n');
    
    if ~isempty(csoc)                                                       % if we have actually performed social moves
        
        % splitting feasible and non feasible solutions
        ffeas = fsoc(csoc<=0,:);
        xfeas = xsoc(csoc<=0,:);
        cfeas = csoc(csoc<=0,:);
        
        fnfeas = fsoc(csoc>0,:);
        xnfeas = xsoc(csoc>0,:);
        cnfeas = csoc(csoc>0,:);
        
        if ~isempty(ffeas)
            
            domtmp=dominance(ffeas,0);                                          % compute dominance on ffeas
            x_tmp=xfeas(domtmp==0,:);                                           % and resize x_tmp keeping only non-dominated elements of domtmp
            f_tmp=ffeas(domtmp==0,:);                                           % idem with f
            c_tmp=cfeas(domtmp==0,:);                                           % and c
            
            [domA,domB]=dominant(memory(memory(:,end)<=0,lx+1:lx+mfit),f_tmp);                 % compute dominance of memory wrt f_tmp, and viceversa
            
            [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
            
            if sum(domB==0)>0                                                   % if there are non-dominated elements in f wrt to memory
                
                qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
                
                [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
                
            end
            
        end
        
        if ~isempty(fnfeas)
            
            domtmp=dominance(cnfeas,0);                                         % compute dominance on ftrial
            x_tmp=xnfeas(domtmp==0,:);                                          % and resize x_tmp keeping only non-dominated elements of domtmp
            f_tmp=fnfeas(domtmp==0,:);                                          % idem with f
            c_tmp=cnfeas(domtmp==0,:);                                          % and c
            
            [domA,domB]=dominant(memory(:,end),c_tmp);                      % compute dominance of memory wrt cf_tmp, and viceversa
            [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
            
            if sum(domB==0)>0                                                       % if there are non-dominated elements in f wrt to memory
                
                qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
                
                [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
                
            end
            
        end
        
    end
    
    oldz = z;
    z=min(memory(:,lx+1:lx+mfit),[],1);
    zstar=max([memory(:,lx+1:lx+mfit); f],[],1);
    if any(z<oldz)
        
        fprintf('SOCIAL MOVES IMPROVED MINIMA\n');
        
        for i=1:mfit
            
            if z(i)<oldz(i)
                
                fprintf('Min of objective %d has changed! Old value = %g, new value = %g \n',i,oldz(i),z(i));
                
            end
            
        end
        
    end
    
    %% FORCED POSITIONING
    
    %     The actual new social repositioning. At first only the agents
    %     associated to the orthogonal subproblems are moved, then, when the
    %     archive is large enough, more and more agents are moved, up to
    %     n_social or the current size of the archive. The "REPOSITIONING"
    %     algorithm proposed above is similar to this one, but should work
    %     better if the objectives are not scaled.
    %
    %                 avail = memory;
    %
    %                 if size(memory,1)>=mfit && n_social>0
    %
    %                     n_move = mfit;
    %
    %                     if size(memory,1)>mfit
    %
    %                         n_move = min(size(memory,1),n_social);
    %
    %                     end
    %
    %                     for i = 1:n_move %problems
    %
    %                         sel = 0;
    %                         val = g_fun(f(id_pop_act_subpr(i),:),lambda(act_subpr(i),:),z,zstar);
    %
    %                         for j = 1:size(avail,1) %new points
    %
    %                             if (avail(j,end)<=0)                                                    % if this social action is in feasible region
    %
    %                                 if (g_fun(avail(j,lx+1:lx+mfit),lambda(act_subpr(i),:),z,zstar)<=val)%||all(avail(j,lx+1:lx+mfit)<=f(id_pop_act_subpr(i),:)) % if sampled position for agent solving i-th subproblem is good (i.e it's gfun is better than the gfun for the old agent)
    %
    %                                     sel = j;
    %                                     val = g_fun(avail(j,lx+1:lx+mfit),lambda(act_subpr(i),:),z,zstar);
    %
    %                                 end
    %
    %                             elseif (avail(j,end)>0)&&(avail(j,end)<cid(id_pop_act_subpr(i)))            % if this social action is not in feasible region but it improves previous constraint violation
    %
    %                                 sel = j;
    %                                 val = g_fun(avail(j,lx+1:lx+mfit),lambda(act_subpr(i),:),z,zstar);
    %
    %
    %                             end
    %
    %                         end
    %
    %                         if sel>0
    %
    %                             f(id_pop_act_subpr(i),:)=avail(sel,lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
    %                             x(id_pop_act_subpr(i),:)=avail(sel,1:lx);                                 % idem with it's position
    %                             cid(id_pop_act_subpr(i))=avail(sel,end);
    %                             avail(sel,:) = [];
    %
    %                             patdirs(id_pop_act_subpr(i)).avail = 1:lx;
    %
    %                             %                 rho(i,1)=options.rhoini;
    %                             %                 rho(i,2)=0;
    %
    %                         end
    %
    %                     end
    %
    %                 end
    
    %% GRADIENT BASED LOCAL MOVES AS FINISHING STEP, ONLY FOR OPTIMAL CONTROL
    
    if refine
        
        if  mod(iter,options.refine_freq)==0 || nfeval>=options.maxnfeval %|| (iter==1)
            
            local_only = 1;
            
            %flip direction for orth subproblems
            int_loc_opt = int_loc_opt+1;
            
            [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho,patdirs,MBH_positions,MADS_dirs]=explore2(memory,x,v,f,cid,nfeval,lambda,act_subpr,id_pop_act_subpr,z,zstar,rho,patdirs,pigr,MBH_positions,MADS_dirs,local_only,int_loc_opt,explore_params);
            
            oldz = z;
            z=min([z;ftrial;discarded.f],[],1);            % update the min of each function (i.e, the columnwise min).
            zstar=max([memory(:,lx+1:lx+mfit);ftrial],[],1);
            
            if any(z<oldz)
                
                fprintf('GRADIENT MOVES IMPROVED MINIMA\n');
                
                for i=1:mfit
                    
                    if z(i)<oldz(i)
                        
                        fprintf('Min of objective %d has changed! Old value = %g, new value = %g \n',i,oldz(i),z(i));
                        
                    end
                    
                end
                
            end
            
            for i=1:options.popsize                                             % for all agents performing local search
                
                if (maxC(i)<=0)                                                     % if current agent is in feasible region
                    
                    if (all(ftrial(i,:)<=f(i,:)))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                        
                        v(i,:) = vtrial(i,:);                                       % given velocity
                        
                        % nornalisation of v, to avoid friction
                        if norm(v(i,:))>0
                            
                            v(i,:) = v(i,:)/norm(v(i,:));
                            
                        end
                        
                        x(i,:)=xtrial(i,:);                                         % position
                        f(i,:)=ftrial(i,:);
                        
                        % objective function value
                        cid(i)=maxC(i);                                             % and max constraint violation
                        patdirs(i).avail = 1:lx;
                        
                    end
                    
                elseif (maxC(i)>0)&&(maxC(i)<cid(i))                                % if current agent is not in feasible region but it's position improves previous constraint violation
                    
                    v(i,:)=xtrial(i,:)-x(i);                                             % store it's velocity
                    x(i,:)=xtrial(i,:);                                             % position
                    f(i,:)=ftrial(i,:);                                             % objective function values
                    cid(i)=maxC(i);                                                 % and max constraint violation
                    
                end
                
            end
            
            % splitting feasible and non feasible solutions
            ffeas = discarded.f(discarded.c<=0,:);
            xfeas = discarded.x(discarded.c<=0,:);
            cfeas = discarded.c(discarded.c<=0,:);
            
            fnfeas = discarded.f(discarded.c>0,:);
            xnfeas = discarded.x(discarded.c>0,:);
            cnfeas = discarded.c(discarded.c>0,:);
            
            if ~isempty(ffeas)
                
                domtmp=dominance(ffeas,0);                                          % compute dominance on ffeas
                x_tmp=xfeas(domtmp==0,:);                                           % and resize x_tmp keeping only non-dominated elements of domtmp
                f_tmp=ffeas(domtmp==0,:);                                           % idem with f
                c_tmp=cfeas(domtmp==0,:);                                           % and c
                
                [domA,domB]=dominant(memory(memory(:,end)<=0,lx+1:lx+mfit),f_tmp);  % compute dominance of memory wrt f_tmp, and viceversa
                
                [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
                
                if sum(domB==0)>0                                                   % if there are non-dominated elements in f wrt to memory
                    
                    qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
                    
                    [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
                    
                end
                
            end
            
            if ~isempty(fnfeas)
                
                domtmp=dominance(cnfeas,0);                                         % compute dominance on ftrial
                x_tmp=xnfeas(domtmp==0,:);                                          % and resize x_tmp keeping only non-dominated elements of domtmp
                f_tmp=fnfeas(domtmp==0,:);                                          % idem with f
                c_tmp=cnfeas(domtmp==0,:);                                          % and c
                
                [domA,domB]=dominant(memory(:,end),c_tmp);                               % compute dominance of memory wrt cf_tmp, and viceversa
                [memory,dd,energy,ener2]=arch_rem(memory,dd,memory(domA~=0,:));
                
                if sum(domB==0)>0                                                       % if there are non-dominated elements in f wrt to memory
                    
                    qq = [x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
                    
                    [memory,dd,energy,ener2,mins,maxs]=arch_shrk6(memory,dd,qq,energy,ener2,mins,maxs,lx,mfit,options.max_arch);
                    
                end
                
            end
            
        end
        
    end
    
    %% REDISTRIBUTE
    
    % Alternative to social move (NOT to the generation of sample but to
    % the actual repositioning of the agent), can be useful especially if
    % objectives are badly scaled
    
    id_agents_orth_subprs = id_pop_act_subpr(1:mfit);   %pick the id of the agents that are currently associated to the orthogonal (first mfit) subproblems
    id_agents_non_orth_subprs = id_pop_act_subpr(mfit+1:end);
    
    avail_memory = memory;
    
    if size(avail_memory,1)>=mfit          % if there are at least mfit points in the archive, satisfy as best as possible the orthogonal subproblems
        
        n_pick = min(size(avail_memory,1),n_social);
        
        % if the archive contains more elements than n_social, take the
        % most well spread from the archive to do the association
        
        if n_pick == n_social
            
            %mm0 = setdiff(memory,avail_memory,'rows','stable');
            
            [mm,ddt,eett,ene2t] = arch_shrk6([],[],memory(1:n_pick,:),0,[],mins,maxs,lx,mfit,n_pick); % first pass, add extremas
            [mm,~,~,~] = arch_shrk6(mm,ddt,memory(n_pick+1:end,:),eett,ene2t,mins,maxs,lx,mfit,n_pick);
            
            avail_memory = mm;%setdiff(mm,mm0,'rows','stable');
            
        end
        
        gtmp = hp_g_fun(avail_memory(:,lx+1:lx+mfit),lambda(act_subpr,:),z,zstar);
        [id_point,id_agent,vals] = associate_agents_positions (gtmp,id_pop_act_subpr);
        
        for i = 1:length(id_point)
            
            % if (all(avail_memory(id_point(i),lx+1:lx+mfit)<=f(id_agent(i),:)))||(any(id_agent(i)==id_pop_act_subpr)&&vals(i,:)<g_fun(f(id_agent(i),:),lambda(act_subpr(id_pop_act_subpr==id_agent(i)),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
            
            f(id_agent(i),:)=avail_memory(id_point(i),lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
            x(id_agent(i),:)=avail_memory(id_point(i),1:lx);                                 % idem with it's position
            cid(id_agent(i))=avail_memory(id_point(i),end);
            
            % else
            
            %     keyboard
            %     fprintf('nope\n')
            % end
            
        end
        
    else
        
        % less than mfit, take the closest in criteria space AMONG the
        % agents solving one individual problem
        
        targets = avail_memory(:,lx+1:lx+mfit);
        candidates = f(id_agents_orth_subprs,:);
        xcandidates = x(id_agents_orth_subprs,:);
        
        % need to check non overlapping agents (if a target is occupied by
        % a non orth agent, an orth agent will be displaced on top of the
        % non orth)
        
        others = f(id_agents_non_orth_subprs,:);
        xothers = x(id_agents_non_orth_subprs,:);
        %         [overlaps,it,io] = intersect(targets,others,'rows','stable');
        %         xoverlaps = x
        %
        %         if any(overlaps)
        %
        %             keyboard;
        %             selection = get_closest(overlaps,candidates);
        %
        %             % change originals
        %
        %             f(id_agents_orth_subprs(selection(:,2)),:) = overlaps(selection(:,1),:);
        %             x(id_agents_orth_subprs(selection(:,2)),:) = xothers(selection(:,1),:);
        %
        %             f(id_agents_non_orth_subprs(selection(:,1)),:) = overlaps(selection(:,2),:);
        %             x(id_agents_non_orth_subprs(selection(:,1)),:) = xothers(selection(:,2),:);
        %
        %             % update candidates
        %             candidates = f(id_agents_orth_subprs,:);
        %
        %         end
        
        selection = get_closest(targets,candidates);
        
        f(id_agents_orth_subprs(selection(:,2)),:)=avail_memory(selection(:,1),lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
        x(id_agents_orth_subprs(selection(:,2)),:)=avail_memory(selection(:,1),1:lx);                                 % idem with it's position
        cid(id_agents_orth_subprs(selection(:,2)))=avail_memory(selection(:,1),end);
        
    end
    
    %% REDISTRIBUTE AND ADAPT SUBPROBLEMS
    
    % Alternative to social move (NOT to the generation of sample but to
    % the actual repositioning of the agent), can be useful especially if
    % objectives are badly scale. Following implementation is buggy, but
    % the idea is promising
    
    %     id_agents_orth_subprs = id_pop_act_subpr(1:mfit);   %pick the id of the agents that are currently associated to the orthogonal (first mfit) subproblems
    %
    %     avail_memory = memory;
    %rng(seed);
    
    %     if size(avail_memory,1)>=mfit          % if there are at least mfit points in the archive, satisfy as best as possible the orthogonal subproblems
    %
    %         % deal with orthogonal subproblems first
    %
    %         for i = id_agents_orth_subprs %for each orth subproblem, move the agent that is solving it to the best position possible
    %
    %             gtmp = hp_g_fun(avail_memory(:,lx+1:lx+mfit),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)';
    %             [ming,posming] = min(gtmp);
    %
    %             if (ming<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) || (all(avail_memory(posming,lx+1:lx+mfit)<=f(i,:)) && (norm(avail_memory(posming,1:lx)-x(i,:))>0) )
    %
    %                 if norm(f(i,:)-avail_memory(posming,lx+1:lx+mfit))>0
    %
    %                     %v(i,:)=0*v(i,:);
    %                     f(i,:)=avail_memory(posming,lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
    %                     x(i,:)=avail_memory(posming,1:lx);                                 % idem with it's position
    %                     cid(i)=avail_memory(posming,end);
    %                     %rho(i,1) = options.rhoini;
    %                     %rho(i,2) = 0;
    %
    %                 end
    %
    %             end
    %
    %         end
    %
    %         [~,idava] = setdiff(avail_memory(:,1:lx),x(id_agents_orth_subprs,:),'rows','stable');
    %         avail_memory = avail_memory(idava,:);
    %
    %         if (size(avail_memory,1) == size(memory,1))   %ugly sanity check
    %
    %             idmin = zeros(mfit,1);
    %
    %             for i = 1:mfit
    %
    %                 [~,idmin(i)] = min(memory(:,lx+i));
    %
    %             end
    %
    %             avail_memory(idmin,:) = [];
    %
    %         end
    %
    %         %if size(avail_memory,1)>=(n_social-mfit)  % more memories than remaining agents, free choice
    %         if size(avail_memory,1)>=mfit  % more memories remaining
    %
    %             n_pick = min(size(avail_memory,1),n_social-mfit);
    %
    %             if n_pick == n_social-mfit
    %
    %                 mm0 = setdiff(memory,avail_memory,'rows','stable');
    %
    %                 [mm,ddt,eett,ene2t] = arch_shrk6([],[],mm0,0,[],mins,maxs,lx,mfit,mfit); % first pass, add extremas
    %                 [mm,~,~,~] = arch_shrk6(mm,ddt,avail_memory,eett,ene2t,mins,maxs,lx,mfit,n_pick+mfit);
    %
    %                 avail_memory = setdiff(mm,mm0,'rows','stable');
    %
    %             end
    %
    %             id_agents_non_orth_subprs = id_pop_act_subpr(mfit+1:end);
    %
    %             gtmp = hp_g_fun(avail_memory(:,lx+1:lx+mfit),lambda(act_subpr(mfit+1:end),:),z,zstar);
    %             [id_point,id_agent,vals] = associate_agents_positions (gtmp,id_agents_non_orth_subprs);
    %
    %             for i = 1:length(id_point)
    %
    %                 oldval = f(id_agent(i),:);
    %
    %                 %if vals(i)< g_fun(oldval,lambda(act_subpr(id_pop_act_subpr==id_agent(i)),:),z,zstar)
    %
    %                     f(id_agent(i),:)=avail_memory(id_point(i),lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
    %                     x(id_agent(i),:)=avail_memory(id_point(i),1:lx);                                 % idem with it's position
    %                     cid(id_agent(i))=avail_memory(id_point(i),end);
    %                     lambda(act_subpr(id_pop_act_subpr==id_agent(i)),:) = (f(id_agent(i),:)-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0))/norm((f(id_agent(i),:)-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0)));
    %
    %                 %end
    %
    %                 %                 if norm(f(id_agent(i),:)-oldval)>0
    %                 %
    %                 %                     v(id_agent(i),:) = 0*v(id_agent(i),:);
    %                 %                     rho(id_agent(i),1) = options.rhoini;
    %                 %                     rho(id_agent(i),2) = 0;
    %                 %
    %                 %                 end
    %
    %             end
    %
    %             %             for i = id_agents_non_orth_subprs
    %             %
    %             %                 gtmp = hp_g_fun(avail_memory(:,lx+1:lx+mfit),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)';
    %             %                 [ming,posming] = min(gtmp);
    %             %
    %             %                 oldval = f(i,:);
    %             %
    %             %                 f(i,:)=avail_memory(posming,lx+1:lx+mfit);                                 % fsamp becomes f(active_problem)
    %             %                 x(i,:)=avail_memory(posming,1:lx);                                 % idem with it's position
    %             %                 cid(i)=avail_memory(posming,end);
    %             %                 avail_memory(posming,:) = [];
    %             %
    %             %                 if norm(f(i,:)-oldval)>0
    %             %
    %             %                     v(i,:)=0*v(i,:);
    %             %                     rho(i,1) = options.rhoini;
    %             %                     rho(i,2) = 0;
    %             %
    %             %                 end
    %             %rng(seed);
    
    %             %             end
    %
    %         end
    %
    %         %else
    %
    %     end
    %
    %     %     if length(unique(f(id_pop_act_subpr,:),'rows','stable'))<n_social
    %     %
    %     %         keyboard
    %     %
    %     %     end
    %
    
    %% SUBPROBLEMS UPDATErng(seed);
    
    
    if n_social&&(mod(iter,supbr_upd_freq)==0)&&all(~isinf(z))                               % if there are active subproblems and they've been active for the last subpr_upd_freq iterations
        
        [act_subpr,id_pop_act_subpr,pigr,lambda_f_best]=upd_act_subpr(pigr,lambda,lambda_f_best,x,f,memory,z,n_social,zstar); % update subproblems
        id_local=setdiff(1:options.popsize,id_pop_act_subpr);               % set the IDs of agents which are performing local actions (i.e differentiate them from the agents which will tacle subproblems through social actions)
        
    end
    
    %% TEMPORARY (FASTER) PLOTTING OF FRONT
    
    if max(memory(:,end)<=0)
        if mfit==3
            plot3(memory(:,lx+1),memory(:,lx+2),memory(:,lx+3),'b.',f(:,1),f(:,2),f(:,3),'r+')
            view(-160,30)
        else
            if mfit==2
                %subplot(2,2,1)
                %plot(discarded.f(:,1),discarded.f(:,2),'k*',fsoc(:,1),fsoc(:,2),'m*')
                plot(memory(:,lx+1),memory(:,lx+2),'b.',f(:,1),f(:,2),'r+',fold(:,1),fold(:,2),'go')
                %hold on
                %line([fold(:,1)'; f(:,1)'],[fold(:,2)'; f(:,2)'],'Color','b')
                %hold off
                %subplot(2,2,2)
                %plot(memory(:,1),memory(:,2),'b.',x(:,1),x(:,2),'r+')%,xtrial(:,1),xtrial(:,2),'g+')
                %subplot(2,2,3)
                %plot(memory(:,3),memory(:,4),'b.',x(:,3),x(:,4),'r+')%,xtrial(:,3),xtrial(:,4),'g+')
                %subplot(2,2,4)
                %plot(memory(:,5),memory(:,6),'b.',x(:,5),x(:,6),'r+')%,xtrial(:,5),xtrial(:,6),'g+')
                %             else
                %
                %                 subplot(1,2,1)
                %                 plot(1:options.popsize,f(:,1),'r+',1:options.popsize,ones(1,options.popsize)*memory(1,lx+1),'b');
                %                 subplot(1,2,2)
                %                 plot(history(:,1),history(:,end-2),'b')
            end
        end
        drawnow
    end
    
    %% ARCHIVE SHRINKING
    
    % NOT NEEDED ANYMORE. With current archiving strategy, the archive
    % always has at most max_arch elements. Could be reactivated, if we
    % really want the inner archive to be larger than the final one. In
    % that cases, this shinking has to be performed in a slightly different
    % way
    
    %      if (archsize>tmp_max_arch)%&&first_shrink                                % if archive is bigger than max allowed within loop, and it's requested to performi a first time shrinking
    %
    %          [memory,qq]=arch_shrk3(memory,lx,mfit,options.max_arch);                      % srhink memory
    %          ener = [ener;qq/size(memory,1)];
    %          archsize=options.max_arch;                                              % update size
    %          %first_shrink=0;                                                     % and never do this in-loop shrinking again
    %          subplot(2,1,1)
    %          plot(memory(:,lx+1),memory(:,lx+2),'bo')
    %          axis equal
    %          subplot(2,1,2)
    %          semilogy(ener)
    %          drawnow
    %      end
    
    %% OUTPUT ON SCREEN
    
    if options.draw_flag>=1
        
        subplot(2,2,2)
        plot(0:1,[options.p_social options.p_social])
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
                set(AX(2),'YLim',[0 options.rhoini],'YGrid','on')
                axis(AX(2),[0 options.popsize 0 0.5])
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
                if size(id_pop_act_subpr)>0
                    
                    plot(f(id_pop_act_subpr(1:2),1),f(id_pop_act_subpr(1:2),2),'gd')
                    
                    if size(id_pop_act_subpr)>2
                        
                        plot(f(id_pop_act_subpr(3:end),1),f(id_pop_act_subpr(3:end),2),'r*')
                        
                    end
                    
                end
                
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
        fprintf('\tMax Constraint violation : %6.3f\n', max(memory(:,end)));
    end
    
    eltime = toc;
    
    if mfit==1
        
        fprintf('Iter: %d, Total nfeval: %d, Elapsed time: %f, fun_evals: %d, time/num_fun_evals: %f\n',iter, nfeval,eltime,nfeval-nfeval_old,eltime/(nfeval-nfeval_old));
        
    else
        
        fprintf('Iter: %d, Arch size: %d, Total nfeval: %d, Elapsed time: %f,fun_evals: %d, time/num_fun_evals: %f, Energy: %f\n',iter, size(memory,1), nfeval,eltime,nfeval-nfeval_old,eltime/(nfeval-nfeval_old),energy);
        
    end
    
    history = [history; iter*ones(size(memory,1),1) memory];
end

%% ARCHIVE SHRINKING BEFORE OUTPUT

if archsize>options.max_arch
    %size(memory,1)
    %     plot3(memory(:,2),memory(:,3),memory(:,4),'bo')
    %     hold on
    [memory,qq]=arch_shrk3(memory,lx,mfit,options.max_arch);
    ener = [ener;qq/size(memory,1)];
    subplot(2,1,1)
    plot(memory(:,lx+1),memory(:,lx+2),'bo')
    axis equal
    subplot(2,1,2)
    semilogy(ener)
    drawnow
    %     plot3(memory(:,2),memory(:,3),memory(:,4),'r*')
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


