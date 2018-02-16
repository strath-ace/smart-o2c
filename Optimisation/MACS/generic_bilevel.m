function [fval,cval,xcorr] = generic_bilevel(x,fhandle,chandle,options)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generic bilevel formulation 

xcorr = x;
cval = [];

x_norm = NORMALIZZA INPUT;

if ~isempty(c)

    [xcorr_norm,~,~,output] = fmincon(@(x) feas_only2(x,x_in_norm),x_norm,[],[],[],[],norm_lb,norm_ub,@(x) chandle,options.fminconoptions);
    
end


% clip vars

x_sol(x_sol<nlb') = nlb(x_sol<nlb');
x_sol(x_sol>nub') = nub(x_sol>nub');

if output.constrviolation<fminconoptions.TolCon
    
    fval =  fhandle(x_sol,problem,0);
    xcorr = xcorr_norm.*problem.scales.scale_opt'; % denormalise x_sol, MACS stores non normalised values in the archive

else
    
    fval = [problem.structure{1}.tf_bounds(2)+1e7*output.constrviolation problem.structure{1}.xf_bounds(2,2)+1e7*output.constrviolation];
    xcorr = xcorr_norm.*problem.scales.scale_opt'; % denormalise x_sol, MACS stores non normalised values in the archive
    
end

end