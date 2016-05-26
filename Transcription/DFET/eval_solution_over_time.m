function [xt,ut] = eval_solution_over_time(x,u,t_0,t_f,t_plot,structure)

ntimes = length(t_plot);
xt = zeros(ntimes,structure.num_eqs);
ut = zeros(ntimes,structure.num_controls);

real_els = t_0+structure.els*(t_f-t_0);

for i = 1:ntimes
    
    % locate element containing current time
    id_elem = 1:structure.num_elems;
    id_elem = min(id_elem((real_els(:,1)<=t_plot(i)).*(real_els(:,2)>=t_plot(i))==1));   
    
    % rescale position of current time in the [-1:1] range for current
    % element
    
    tmax = real_els(id_elem,2);
    tmin = real_els(id_elem,1);
    
    norm_time = 2*(t_plot(i)-(tmax+tmin)/2)/(tmax-tmin);
    
    xt(i,:) = (structure.state_basis{id_elem}(norm_time)*x(:,id_elem))';
    ut(i,:) = (structure.control_basis{id_elem}(norm_time)*u(:,id_elem))';

end

end