function [value,isterminal,direction] = BallodeEventsFcn(t,y)
% Locate the time when height passes through zero in a
% decreasing direction and stop integration.
value      = [y(1)]; % Detect height = 0 or 0.5
isterminal = [ 1];         % Stop the integration
direction  = [-1];         % Down 