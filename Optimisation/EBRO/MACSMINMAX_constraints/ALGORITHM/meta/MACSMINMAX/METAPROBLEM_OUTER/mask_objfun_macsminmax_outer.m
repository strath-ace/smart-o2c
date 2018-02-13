function [masked] = mask_objfun_macsminmax_outer(d,par_objfun)

% for i=1:length(d)
%     if d(i)>1
%         a=1;
%     end
% end
[masked, ~, ~] = u_validation( par_objfun.problem_fix_d, d, par_objfun.u_record, par_objfun.local_search, par_objfun.objectives);

return
