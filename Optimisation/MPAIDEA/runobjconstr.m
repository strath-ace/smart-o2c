function [x,f,eflag,outpt] = runobjconstr(x0, fc, flag_LG, opts, LB, UB, varargin)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% runobjconstr
% This function is used to:
% - compute objective function including weighted penalty for constrained
% problem to be handled by the DE in MPAIDEA
% - compute objective function for unconstrained problems for DE of MPAIDEA
% - compute solution of local search using fmincon, for both constrained
% and unconstrained problem
% Refer to: https://uk.mathworks.com/help/optim/ug/objective-and-nonlinear-constraints-in-the-same-function.html
% 
%% Inputs:
% 
% * x0: initial guess for fmincon, or individual of the population where to
%       evaluate the objective function
% * fc: structure containing information about the objective and
%       constraints, namely:
%       * fc.obj: objective function handle 
%       * fc.constr: constraints function handle
%       * fc.w_ceq: weights for equality constraints
%       * fc.w_c: weights for inequality constraints
%       * fc.obj_constr: flag used to specify if the objectives and the 
%                        constraints are defined in the same function 
%                        (fc.obj_constr = 1) or in different functions 
%                        (fc.obj_constr = 0)
%       * flag_LG: structure defining if this function should be used for 
%                  the local search of MP-AIDEA (flag_LG.local = 1, 
%                  flag_LG.global = 1) or for the DE evaluation in MP-AIDEA
%                  (flag_LG.local = 0, flag_LG.global = 1)
% * opts: options for optimisation with fmincon
% * LB   -> lower boundaries for fmincon
% * UB   -> upper boundaries for fmincon
% * varargin -> additional inputs

%% Outputs:
% 
% * x: solution vector of fmincon if flag_LG.local=1 and flag_LG.global =
%      0; empty if flag_LG.local = 0 and flag_LG.global = 1
% * f: solution value of fmincon or weighted objective function for DE
% * eflag: exitflag for fmincon if flag_LG.local=1 and flag_LG.global = 0
% * outpt: output of fmincon if flag_LG.local=1 and flag_LG.global = 0

%% Author(s): Marilena Di Carlo (2015)
% email: marilena.di-carlo@strath.ac.uk


if nargin == 1 % No options supplied
    opts = [];
end

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


% -------------------------------------------------------------------------
% Local search: fmincon
% -------------------------------------------------------------------------
if flag_LG.local == 1
    
    fun = @objfun; % the objective function, nested below
%     fun = @(xs)fc.obj(x,varargin);
    
    % Problem with upper and lower boundaries
    if ~isempty(LB) && ~isempty(UB)
        
        % If problem has constraints:
        if ~isempty(fc.constr)
            cfun = @constr; % the constraint function, nested below
            % Call fmincon
            [x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],LB, UB,cfun,opts,varargin{:});
        else
            % Call fmincon
            [x,f,eflag,outpt] = fmincon(fc.obj,x0,[],[],[],[],LB, UB,[],opts,varargin{:});
        end
        
    % Unbounded problem:
    else
        % Call fmincon
        [x,f,eflag,outpt] = fminunc(fun,x0,opts);
    end

% -------------------------------------------------------------------------
% Global search DE: ojective function with penalty
% -------------------------------------------------------------------------
elseif flag_LG.global == 1
    
    x = x0;
    
    % Evaluate objective function
    yy = objfun(x, varargin{:});
    
    % If problem is constrained, evaluate function of constraints
    if ~isempty(fc.constr)
        [c,ceq] = constr(x, varargin{:});
    else
        c = [];
        ceq = [];
    end
    
    if (fc.weighted && (~isempty(ceq) || ~isempty(c)) ) || ...
            ( isempty(ceq) && isempty(c) )
        
        % Defined weighted objective function:
        % Both equality and inequality constraints
        if ~isempty(ceq) && ~isempty(c)
            f = yy + fc.w_ceq * norm(ceq)  + fc.w_c * ( abs(c) .* (c>0 == 1));
            % Only equality constraints
        elseif ~isempty(ceq) && isempty(c)
            f = yy + fc.w_ceq * norm(ceq)  ;
            % Only inequality constraints
        elseif isempty(ceq) && ~isempty(c)
            f = yy + fc.w_c * ( abs(c) .* (c>0 == 1));
            % No constraints
        elseif isempty(ceq) && isempty(c)
            f = yy ;
        end
        
    elseif fc.weighted == 0 && (~isempty(ceq) || ~isempty(c))
        
        f.yy = yy;
        f.ceq = norm(ceq);
        f.c = c;
        
        if norm(f.ceq) > fc.ceq_eps  && any(f.c > 0)
            f.non_feas_ceq = 1;
            f.non_feas_c    = 1;
        elseif norm(f.ceq) > fc.ceq_eps  && all(f.c < 0)
            f.non_feas_ceq = 1;
            f.non_feas_c    = 0;
        elseif norm(f.ceq) < fc.ceq_eps  && any(f.c > 0)
            f.non_feas_ceq = 0;
            f.non_feas_c    = 1;
        else
            f.non_feas_ceq = 0;
            f.non_feas_c = 0;
        end
        
    end
    
    eflag = [];
    outpt = [];
    
end

    function y = objfun(x, varargin)
        
        % Check if computation is necessary (only if
        % objective function and constraint functions are the same)
        if ~isequal(x,xLast) && fc.obj_constr

            [myf,myc,myceq] = feval(fc.obj,x,varargin{:});
            xLast = x;
            
        elseif isequal(x,xLast) && fc.obj_constr
            
        elseif ~isequal(x,xLast) && ~fc.obj_constr
            % Evaluate objective function            
            [myf] = feval(fc.obj,x,varargin{:});
        end
        
        y = myf;
    end

    function [c,ceq] = constr(x, varargin)
        % Check if computation is necessary
        if ~isequal(x,xLast) && fc.obj_constr 
            
            [myf,myc,myceq] = feval(fc.constr, x,varargin{:});
            xLast = x;
        
        elseif isequal(x,xLast) && fc.obj_constr
            
        elseif ~isequal(x,xLast) && ~fc.obj_constr
            % Evaluate constraint function   
            [myc,myceq] = feval(fc.constr, x,varargin{:});
            
        end
        
        % Now compute constraint functions
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end


end