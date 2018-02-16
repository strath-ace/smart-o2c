function [val,x_sol] = DFET_bilevel(x_in,K,L,maxtrials,problem,fminconoptions)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% normalise input

x_in_norm = x_in'./problem.scales.scale_opt;

% don't allow variations in static parameters, thus remove them from
% solution vector

nlb = problem.norm_lb;
nub = problem.norm_ub;

% if we're sampling on the boundary, move it a bit in the inner side
% because interior point doesn't like initial guesses laying on the bounds
%x_in_norm(x_in_norm==nlb) = nlb(x_in_norm==nlb)+1e-6;
%x_in_norm(x_in_norm==nub) = nub(x_in_norm==nub)-1e-6;

jacflag = 1;

doit = 1;
options = fminconoptions;

x_guess = x_in_norm;

tstart = tic;

[c,ceq] = multiphase_constr_full(problem,x_guess,0,tstart,problem.timeout_feasible);

feas0 = max([abs(ceq);c]);
trials = 0;

if feas0>options.TolCon
    
    while doit  %workaround for stupid behaviour of fmincon (reporting a final feasibility different than the real one)
        
        timeout = 0;
        tstart = tic;
        
        try
            
            [x_sol,~,~,output] = fmincon(@(x) feas_only2(x,x_guess),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,jacflag,tstart,problem.timeout_feasible),options);
            
        catch EX
            
            if strcmp(EX.identifier,'multiphase_constr_full:timeOutReached')
                
                timeout = 1;
                x_sol = x_guess;
                fprintf('Timeout of %f seconds reached\n', problem.timeout_feasible);
                
            else
                
                rethrow(EX);
                
            end
            
        end
        
        [c,ceq] = multiphase_constr_full(problem,x_sol,0,tstart,-1); % with a negative timeout it has all the time it needs. This way we never generate an exception here, where it's not meant to be any time limit
        newfeas = max([abs(ceq);c]);
        
        if ~timeout
            
            if (newfeas<=options.TolCon)
                
                doit = 0;
                
            else
                
                if (output.stepsize>options.StepTolerance) && trials<maxtrials
                    
                    x_guess = x_sol;
                    %feas = newfeas;
                    trials = trials+1;
                    
                else
                    
                    doit=0;
                    
                end
                
            end
            
        else
            
            doit = 0;
            
        end
        
    end
    
else
    
    x_sol = x_guess;
    newfeas = feas0;
    
end

% clip vars

x_sol(x_sol<nlb) = nlb(x_sol<nlb);
x_sol(x_sol>nub) = nub(x_sol>nub);

% if we removed static vars from optimisation, re include them to have the
% full solution vector

if newfeas<options.TolCon
    
    val =  multiphase_objectives (x_sol,problem,0);
    x_sol = x_sol.*problem.scales.scale_opt; % denormalise x_sol, MACS stores non normalised values in the archive
    
else
    
    val = L+K*newfeas;
    x_sol = x_sol.*problem.scales.scale_opt; % denormalise x_sol, MACS stores non normalised values in the archive
    
end


end
