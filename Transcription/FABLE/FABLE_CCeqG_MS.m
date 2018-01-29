function [myf,myc,myceq] =  FABLE_CCeqG_MS(x, parameters, options, constants, inputs)

% =========================================================================
% Function for the simultanous computation of objective function,
% constraints and their gradients for FABLE
% Input: x          -> solution vector
%        parameters -> structure containing the parameters of the problem
%        options    -> structure containing the options for the problem
%        constants  -> structure of constants
% 
% Output: myf   -> structure containing objective function and its gradient
%         myc   -> structure containing inequality constraints and their
%                  gradients
%         myceq -> structure containing equality constraints and their
%                  gradients
% =========================================================================
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk
% =========================================================================

% Compute objective functions and constraints
if strcmp(options.MS_SS,'MS') 
    [J,C,Ceq,GJ, GC, GCeq] = FABLE_transcription_MS(x, parameters, options, constants);
elseif strcmp(options.MS_SS,'SS')
    [J,C,Ceq,GJ, GC,GCeq] = FABLE_transcription_SS(x, parameters, options, constants, inputs);
    myf.GJ = GJ;
end

myf.J = J;
myc.C = C;
myceq.Ceq = Ceq;



% If gradient of constraints or of objective functions is provided by the
% user, compute them with finite differences
if strcmp(options.options_fmincon.GradConstr, 'on') && ...
    strcmp(options.options_fmincon.GradObj, 'on')
    
    myceq.GCeq = GCeq';
    myf.GJ = GJ;
    myc.GC = GC';
    
elseif strcmp(options.options_fmincon.GradConstr, 'on') && ...
    strcmp(options.options_fmincon.GradObj, 'off')
    
    myceq.GCeq = GCeq';
    myc.GC = GC';
    myf.GJ = [];
else
    myf.GJ = [];
    myc.GC = [];
    myceq.GCeq = [];
end