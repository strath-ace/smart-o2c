function [memories, varargout] = optimise(fun,input,options,varargin)

xtrial = (input.vub-input.vlb)/2;
num_objs = length(fun(xtrial,varargin));

oldpath = addpath('MPAIDEA','MACS');

if num_objs <2
        
    % call AMPIDEA
    
    [memories,B_mean,bubble, archivebest,options] = MP_AIDEA_ALR(fun, input, options, varargin);
    
    varargout{1} = B_mean;
    varargout{2} = bubble;
    varargout{3} = archivebest;
    varargout{4} = options;
        
else
    
    % call MACS
   
    [memories,nfeval,ener]=macs7v16OC(fun,[],input.vlb,input.vub,options,[],[],varargin);
    varargout{1} = nfeval;
    varargout{2} = ener;
    
end

rmpath('MPAIDEA','MACS');

end