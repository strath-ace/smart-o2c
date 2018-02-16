function [] = animate_front(mem)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates an animation of the evolution of the Pareto front 

nit = max(mem.history(:,1));

mino1 = min(mem.history(:,end-3));
maxo1 = max(mem.history(:,end-3));
mino2 = min(mem.history(:,end-2));
maxo2 = max(mem.history(:,end-2));
deltao1 = maxo1-mino1;
deltao2 = maxo2-mino2;

figure()

for i = 1:nit
   
    plot(mem.history(mem.history(:,1)==i,end-3),mem.history(mem.history(:,1)==i,end-2),'b*')
    axis([mino1-deltao1*0.2 maxo1+deltao1*0.2 mino2-deltao2*0.2 maxo2+deltao2*0.2])
    pause(0.5)
    drawnow
    
end
