function [x, f, eflag, outpt] = FABLE_runobjconstr(x0,parameters,options,constants,inputs,A,b,Aeq, beq, LB,UB)

% =========================================================================
% Function for the call to fmincon for FABLE.
% Used for cases when objective and constraints are defined in the same function. 
% Fmincon would call the same function twice to evaluate them. Using this
% method the function is called instead only once.
%
% From Matlab help: Objective and Nonlinear Constraints in the Same Function
% https://uk.mathworks.com/help/optim/ug/objective-and-nonlinear-constraints-in-the-same-function.html
%  ========================================================================
% Input: x0 - > initial guess for NLP
%        parameters -> structure containing the parameters of the problem
%        options    -> structure containing the options for the problem
%        constants  -> structure of constants
%        A          -> matrix for the linear constraints A*x<=b
%        b          -> array for linear constraints A*x<=b
%        Aeq        -> matrix for linear constraints Aeq*x=beq
%        beq        ->
%        LB         -> lower boundaries of the search space
%        UB         -> upper boundaries of the search space
%
% Output: x         -> optimal solution vector
%         f         -> optimal objective function
%         eflag     -> exitflag of fmincon
%         outpt     -> output of fmincon
%
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%%

xLast = []; % Last place computeall was called
myf = [];   % Use for objective at xLast
myc = [];   % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


fun = @objfun;  % the objective function, nested below
cfun = @constr; % the constraint function, nested below

% Scale variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 = (x0 - LB) ./ (UB - LB);
%
% LB2 = zeros(1,numel(x0));
% UB2 = ones(1,numel(x0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call fmincon
[x,f,eflag,outpt] = fmincon(fun,x0,A, b, Aeq, beq, LB, UB,cfun,options.options_fmincon);


%% Objective function
    function [y,myGf] = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] =  FABLE_CCeqG_MS(x, parameters, options, constants, inputs);
            xLast = x;
        end
        
        y = myf.J;
        myGf = myf.GJ;
        
    end

%% Constraints function
    function [c,ceq,myGc,myGCeq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] =  FABLE_CCeqG_MS(x, parameters, options, constants, inputs);
            
            xLast = x;
        end
        
        % Now compute constraint functions
        c = myc.C; 
        ceq = myceq.Ceq;
        max(abs(ceq));
             
        myGc = myc.GC;
        myGCeq = myceq.GCeq;
               
    end

end