function  Hs = vw2hs(Vw)
% This function converts the average wind speed Vw to significan wave 
% heihgt Hs accodding to: Hs =  0.21 Vw^2 / g 
%
% Inputs:  
%   Vw : Wind speed (m/s)
%
% Outputs:
%   Hs : Significant wave height (m)
%
% References:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
% Author:     Tristan Perez 
% Date:       2005-03-12
% Revisions: 

g = 9.81;
Hs =  0.21 * Vw^2 / g; 
