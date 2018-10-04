function  Hs = vw2hs(Vw)
% This function converts the average wind
% speed Vw  to significan wave heihgt Hs accodding to 
% 
%      Hs =  0.21 Vw^2 / g 
%
% Use:  [Hs]=vw2hs(Vw)
% Vw -Wind speed [m/s]
% Hs -Significant wave height [m].
%
% Reference:
% T.I. Fossen (2002) "Marine Control Systems" Marine Cybernetics. 
%
% Created by: 2005-03-12 Tristan Perez 

g = 9.81;
Hs =  0.21*Vw^2 / g; 