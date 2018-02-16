function [J]=mexspaceartCassiniNOdsmA_MO(ya)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------


ll=[-1000 30   100  30   400   1000
     0    400  470  400  2000  6000];
 
y=(ll(2,:)-ll(1,:)).*ya+ll(1,:);

%clamp!!!!
y(y<ll(1,:)) = ll(y<ll(1,:));
y(y>ll(2,:)) = ll(y>ll(2,:));

J=[mexspaceartCassiniNOdsm(y) sum(y(2:end))];
% 
% if J(1)>15
%     
%     J = [J(1)^2 J(1)^2+5e3];
%     
% end