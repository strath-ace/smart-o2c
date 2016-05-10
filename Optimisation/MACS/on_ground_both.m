function [value,isterminal,direction] = on_ground_both(t,y,erad,toll)

% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.

% if (norm(y(1:3))-erad)<0 || (norm(y(4:6))-erad)<0
%    keyboard 
% end

value = [norm(y(1:3))-erad; norm(y(4:6))-erad ; norm(y(1:3)-y(4:6))-toll];     % Detect height = 0
isterminal = [1;1;1];   % Stop the integration
direction = [0;0;0];   % Negative direction only
%direction = [-1;-1;-1];   % Negative direction only

end