function [value,isterminal,direction] = on_ground(t,y,erad)

% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.

value = norm(y(1:3))-erad;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only

end