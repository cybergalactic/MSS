function  Vw = hs2vw(Hs)
% This function converts the significan wave heihgt Hs into an equivalent wind
% speed Vw accodding to 
% 
%        Vwc= sqrt((Hs*g/0.21));
%
% Use:  [Vw]=hs2vw(Hs)
% Hs -Significant wave height [m].
% Vw -Wind speed [m/s]

g =9.81;
Vw = sqrt((Hs*g/0.21));