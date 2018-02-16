function [x_guess_probl] = generate_guess(problem,adim,varargin)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function generates an initial guess for the whole multiphase
% problem. The phases are initialised sequentially and independently, and
% then are simply appended

% Lorenzo A. Ricciardi, 2017

if adim<0 || adim>1 || mod(adim,1)>0
    
    error('adim must be either 0 or 1');
    
end

if length(varargin)>1
    
    error('Too many input parameters in generate_guess')
    
end

if (~isempty(varargin)) && (~isequal(length(varargin{1}),length(problem.lb)))
    
    error('Supplied guess to generate_guess has wrong dimensions');
    
end

x_guess_probl = [];

%phase_mask = [];

tt = zeros(problem.num_phases+1,1);

for i = 1:problem.num_phases
    
    if isempty(varargin)    % construct guess using the values supplied in this phase definition
        
        x_0 = problem.structure{i}.x_0;
        t_0 = problem.structure{i}.t_0;
        t_f = problem.structure{i}.t_f;
        static = problem.structure{i}.other_vars_guess;
        u_nodes = repmat(reshape(ones((problem.structure{i}.control_order+1),1)*problem.structure{i}.control_guess',problem.structure{i}.num_controls*(problem.structure{i}.control_order+1),1),problem.structure{i}.num_elems,1);
        
    else                    % construct guess using externally generated values (typically from MACS)
        
        this_phase_guess = varargin{1};
        
        if problem.structure{i}.imposed_t0
            
            t_0 = problem.structure{i}.t_0;
            
        else
            
            t_0 = this_phase_guess(problem.phase_mask==i.*problem.t0_vars);
            
        end
        
        if problem.structure{i}.imposed_tf
            
            t_f = problem.structure{i}.t_f;
            
            
        else
            
            t_f = this_phase_guess(problem.phase_mask==i.*problem.tf_vars);
            
        end
        
        x_0 = problem.structure{i}.x_0;
        
        x_0(~problem.structure{i}.imposed_initial_states) = this_phase_guess(problem.phase_mask==i.*problem.x0_vars);
        
        static = this_phase_guess(problem.phase_mask==i.*problem.other_vars);
        static = static(:); % ensure column vector
        u_nodes = this_phase_guess(problem.phase_mask==i.*problem.control_vars);
        
    end
    
    fprintf('Generating feasible guess for phase %d... ',i);
    
    tic
    x_guess = make_first_guess(x_0,t_0,t_f,static,u_nodes,problem.structure{i});
    tt(i) = toc;
    fprintf('%f s\n',tt(i));
    % normalise all
    
    if strcmp(problem.structure{i}.init_type,'CloneIC')
        
        x_guess = [problem.structure{i}.other_vars_guess;x_guess];
        
    end
    x_guess = x_guess./problem.structure{i}.scales.scale_opt;
    
    % hard clipping of x_guess within limits
    
    %x_guess(x_guess<problem.structure{i}.norm_lbv) = problem.structure{i}.norm_lbv(x_guess<problem.structure{i}.norm_lbv);
    %x_guess(x_guess>problem.structure{i}.norm_ubv) = problem.structure{i}.norm_ubv(x_guess>problem.structure{i}.norm_ubv);
        
    %phase_mask = [phase_mask; i*ones(size(x_guess))];
    
    x_guess_probl = [x_guess_probl; x_guess];
    
end

%% Link all phases (BEWARE, TIME ORDERING OF PHASES CAN BE RESHUFFLED BY ABOVE METHOD)

options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',length(problem.norm_lb),'MaxIter',length(problem.norm_lb),'GradConstr','on','GradObj','on','TolFun',1e3,'Tolcon',1e-6,'TolX',1e-12,'DerivativeCheck','off');%,'PlotFcns',@(x,optimValues,state) myoptimplot(x,optimValues,state,structure,t_0,x_0,x_f));

doit = 1;
maxit = 3;
it = 1;
fprintf('Linking phases into fully feasible trajectory... ');
nfunc = 0;
nit = 0;

while doit && (it<maxit)

    tstart = tic;
    timeout = 0;
    fprintf('trial %d, ', it);
    
    try
    
        [x_guess_probl,~,~,output] = fmincon(@(x) feas_only(x),x_guess_probl,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1,tstart,-1),options);
    
    catch EX
        
        if strcmp(EX.identifier,'multiphase_constr_full:timeOutReached')
       
            timeout = 1;
            fprintf('Timeout of %f seconds reached\n', problem.timeout_feasible);
        
        else
            
            rethrow(EX);
            
        end
        
    end
    
    if ~timeout
        
    it = it+1;
    
    [c,ceq] = multiphase_constr_full(problem,x_guess_probl,0,tstart,-1);
    newfeas = max([abs(ceq);c]);
   
    nfunc = nfunc + output.funcCount;
    nit = nit + output.iterations;       
    
    if (newfeas<=options.TolCon) || (output.stepsize<options.StepTolerance)
        
        doit = 0;

    end
    
    else
       
        doit = 0;
        newfeas = inf;
        
    end

end

tt(end) = toc(tstart);

if(newfeas<=options.TolCon)
    
    fprintf('%f s (%d fevals, %d iters)\n',tt(end),nfunc,nit);
    fprintf('Fully feasible solution found in %f s\n\n',sum(tt));

else
    
    fprintf('FAILED!\n\n');
    
end

if ~adim % solution has to be returned in dimensional coordinates
    
    x_guess_probl = x_guess_probl.*problem.scales.scale_opt;
    
end

end