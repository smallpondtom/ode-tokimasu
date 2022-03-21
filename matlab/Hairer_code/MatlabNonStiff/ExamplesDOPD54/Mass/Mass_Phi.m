function [Phi,tLag,tLagTab] = Mass_Phi(t,y,varargin)
Phi          =  [cos(t),0];
tLag         = -2*pi;
tLagTab      = zeros(2,1);
tLagTab(1,1) = tLag;
