function [KT, KQ] = wageningen(Ja,PD,AEAO,z)
% [KT, KQ] = wageningen(Ja,PD,AEAO,z) computes the thrust KT and torque KQ
% coefficients of the Wageningen B-series propellers.
% 
% Inputs:
%   Ja = Va/(n*D)  : Advance number
%   PD             : Pitch/diameter ratio (typically 0.5-2.5)
%   AEAO           : Blade area ratio (ratio of expanded blade area to 
%                  : Propeller disk area)
%   z              : Number of propeller blades
%
% Usage:
%  See exWageningen.m 
%
% The B-series propellers were designed and tested at the Netherlands Ship
% Model Basin in Wageningen. The open_water characteristics of 120
% propeller models of the B_series were tested at N.S.M.B. and analyzed
% with multiple polynomial regression analysis. The derived polynomials
% express the thrust and torque coefficients in terms of the number of
% blades, the blade area ratio, the pitch_diameter ratio and the advance
% coefficient.
%
% Reference: 
%   Barnitsas, M.M., Ray, D. and Kinley, P. (1981).
%   KT, KQ and Efficiency Curves for the Wageningen B-Series Propellers
%   http://deepblue.lib.umich.edu/handle/2027.42/3557
%
% Author:    Thor I. Fossen
% Date:      7 October 2018
% Revisions: 
%   18 June 2019 - minor updates

load('WageningData.mat');
    
KT = sum(WagCThrust_stuv.*((Ja).^WagThrust_s).*(PD.^WagThrust_t).*...
    (AEAO.^WagThrust_u).*(z.^WagThrust_v));
KQ = sum(WagCTorque_stuv.*((Ja).^WagTorque_s).*(PD.^WagTorque_t).*...
    (AEAO.^WagTorque_u).*(z.^WagTorque_v));
