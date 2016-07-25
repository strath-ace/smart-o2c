% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [c,ceq,Jc,Jceq] = dynamics(structure,x_in,x_0,x_f,t_0,t_f,jacflag)

c = [];
Jc = [];

[ceq,Jceq] = eval_constraints(structure.f,structure,x_in,x_0,x_f,t_0,t_f,structure.uniform_els,jacflag,structure.dfx,structure.dfu);

Jceq = Jceq';

end