function [value,isterminal,direction] = BallodeEventsFcn(t,y)
% Locate the time when height passes through zero in a
% decreasing direction and stop integration.
% Locate the time when height passes 5 in each direction and stop 
% integration.
value      = [y(1), y(1)-5]; % Detect height = 0 or 0.5
isterminal = [ 1,1];         % Stop the integration
direction  = [-1,0];         % Down for zero, each direction for height=5