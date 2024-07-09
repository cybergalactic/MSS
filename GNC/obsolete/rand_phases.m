function [ph]=rand_phases(nc)
%Function to generate a uniformely 
%distributed vector of random phases 
%in the interval [-pi pi].
%
% Use: [ph]=rand_phases(nc)
%
% ph -vector of random phases 
% nc- number of components
% 
% Created by: Tristan Perez in 2001  
% Last mod. by: Tristan Perez 
% Date: 9 March 2005

ph = 2*pi * (rand([nc 1])-.5);
