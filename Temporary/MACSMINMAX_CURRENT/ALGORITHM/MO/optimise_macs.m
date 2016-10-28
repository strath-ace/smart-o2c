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
% * x : archive containing the solution vectors
% * fval : archive containing obejctive function values
% * exitflag : UNUSED AT MOMENT
% * output : additional outputs from MACS
%          * memory : the entire archive
%          * nfeval : effective number of objective function evaluations
%          * ener : value of the energy of the archive, used by Energy
%                   Based Archiving strategy
%
% Author: Lorenzo A. Ricciardi (2015-2016)
% email: lorenzo.ricciardi@strath.ac.uk

% Run MACS
[memory,nfeval,ener]=macs7v16OC(fitnessfcn,[],LB,UB,options,[],[],varargin{1:end});

% Output
x    = memory(:,1:length(LB));
fval = memory(:,length(LB)+1:end-2);

output.memory               = memory;
output.nfeval               = nfeval;
output.ener                 = ener;
exitflag                    = 1;

end
