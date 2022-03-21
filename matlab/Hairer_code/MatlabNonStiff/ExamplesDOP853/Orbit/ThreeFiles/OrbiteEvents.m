function [value,isterminal,direction] = OrbitEvents(t,y,mu,mustar,y0)
% Locate the time when the object returns closest to the
% initial point y0 and starts to move away; stop integration.
% Also locate the time when the object is farthest from the 
% initial point y0 and starts to move closer.
% 
% The current distance of the body is
% 
%   DSQ = (y(1)-y0(1))^2 + (y(2)-y0(2))^2 
%       = <y(1:2)-y0(1:2),y(1:2)-y0(1:2)>
% 
% A local minimum of DSQ occurs when d/dt DSQ crosses zero 
% heading in the positive direction. Compute d(DSQ)/dt as
% 
%  d(DSQ)/dt = 2*(y(1:2)-y0(1:2))'*dy(1:2)/dt = ... 
%                 2*(y(1:2)-y0(1:2))'*y(3:4)
% 
dDSQdt = 2 * ((y(1:2)-y0(1:2))' * y(3:4));
value = [dDSQdt; dDSQdt];
isterminal = [1; 0];            % Stop at local minimum
direction = [1; -1];            % [local minimum, local maximum]

