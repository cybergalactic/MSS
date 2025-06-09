function Cd_2D = cylinderDrag(L,B,nu_r)
% Reynolds number and aspect ratio dependent 2D drag coefficient for cylinders.
% Literature data was digitized and interpolation is used to compute intermediate points. 
%
%  Cd_2D = cylinderDrag(L,B,nu_r)
%
% Output: 
%    Cd_2D: 2D drag coefficient
%
% Inputs:
%    L: length
%    B: beam 
%    nu_r = [u-u_c, v-v_c, w-w_c, p, q, r]': relative velocity vector
%
% Author: M. Seidl
% Date:   09 Jun 2025
%
% Reference:
% DNV-RP-C205, Environmental conditions and environmental loads
% Recommended practice, Edition 2025-04

% CD_DATA = [Re  Cd]
CD_DATA = [...
10211.0405297256 1.20769
15543.5423211490 1.20369
24434.8681540938 1.20823
35719.8807222993 1.21060
60880.7508493947 1.21093
86174.7436190610 1.20901
118086.431856986 1.21134
149262.255516420 1.21148
214708.876786124 1.21171
230940.445891927 1.20109
271581.385758833 1.16918
297098.185911873 1.11163
325252.018159382 1.00926
361942.848517458 0.89411
409798.156149477 0.70855
475302.392082697 0.53155
524734.960427919 0.40785
578921.928291835 0.32683
638292.879987135 0.28422
744046.366191920 0.29711
853237.628592033 0.32280
1068685.70450627 0.38909
1393381.37138896 0.47033
1831850.38455424 0.53878
2258395.50114555 0.58372
2899011.91287642 0.62868
3663279.03082383 0.64803
4937305.49552761 0.67808];

% kappa is the aspect ratio (AR) dependent drag reduction factor which is 
% different for the sub-critical (Re < 2e5) and super-critical (Re >= 2e5) regime.

% KAPPA_SUBCRITICAL_DATA = [AR  kappa]
KAPPA_SUBCRITICAL_DATA = [...
  2 0.58
  5 0.62
 10 0.68
 20 0.74
 40 0.82
 50 0.87
100 0.98];

% KAPPA_SUPERCRITICAL_DATA = [AR  kappa]
KAPPA_SUPERCRITICAL_DATA = [...
  2 0.80
  5 0.80
 10 0.82
 20 0.90
 40 0.98
 50 0.99
100 1.00];

U_crossflow = sqrt(nu_r(2)^2+nu_r(3)^2); % cross-flow velocity
Re = U_crossflow * L * 1e6; % Reynolds number for water

% Cd interpolation
if Re < CD_DATA(1,1)
    Cd = CD_DATA(1,2);
elseif Re > CD_DATA(end,1)
    Cd = CD_DATA(end,2);
else
    Cd = interp1(CD_DATA(:,1),CD_DATA(:,2),Re);
end

% Aspect ratio dependent reduction
if Re < 2e5
    if L/B < KAPPA_SUBCRITICAL_DATA(1,1)
        kappa = KAPPA_SUBCRITICAL_DATA(1,2);
    elseif L/B > KAPPA_SUBCRITICAL_DATA(end,1)
        kappa = KAPPA_SUBCRITICAL_DATA(end,2);
    else
        kappa = interp1(KAPPA_SUBCRITICAL_DATA(:,1),KAPPA_SUBCRITICAL_DATA(:,2),L/B);
    end
else
    if L/B < KAPPA_SUPERCRITICAL_DATA(1,1)
        kappa = KAPPA_SUPERCRITICAL_DATA(1,2);
    elseif L/B > KAPPA_SUPERCRITICAL_DATA(end,1)
        kappa = KAPPA_SUPERCRITICAL_DATA(end,2);
    else
        kappa = interp1(KAPPA_SUPERCRITICAL_DATA(:,1),KAPPA_SUPERCRITICAL_DATA(:,2),L/B);
    end
end

Cd_2D = Cd * kappa;
