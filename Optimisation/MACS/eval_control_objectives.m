function val = eval_control_objectives(x,structure,x_f,t_0,t_f)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%

[x_out,u,x_b] = extract_solution(x,structure,x_f);       
     
val = eval_cost_functions(structure.g,structure.weights,x_out,u,x_b,[t_0 t_f],structure.uniform_els,structure,0,[],[])';

end
