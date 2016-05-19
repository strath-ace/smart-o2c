function [x,fval,exitflag,output] = optimise_macs(fitnessfcn,LB,UB,options,varargin)

%% optimise_macs

%% Inputs:
%
% * fitnessfcn : function handle to function to optimise
% * LB : lower boundaries
% * UB : upper boundaries
% * options : 
% * varargin : additional parameters on which fitnessfcn depends but that
% are NOT optimisation variables
%
%% Output:
% * y : explanation
%
% Author: Lorenzo A. Ricciardi
% email: lorenzo.ricciardi@strath.ac.uk


% Run MACS
[memory,nfeval,ener]=macs7v16OC(fitnessfcn,[],LB,UB,options,[],[],varargin{1:end});

% Output
x    = memory(:,1:length(LB));
fval = memory(:,length(LB):end-2);


output.memory               = memory;
output.nfeval               = nfeval;
output.ener                 = ener;
exitflag                    = 1;

end