function  Hs = vw2hs(Vw)
% vw2hs is compatibel with MATLAB and GNU Octave (www.octave.org).
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
%   T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
%       Control, 2nd edition, John Wiley & Sons. Ltd., Chichester, UK.
%
% Author:     Tristan Perez 
% Date:       2005-03-12
% Revisions: 

g = 9.81;
Hs =  0.21 * Vw^2 / g; 