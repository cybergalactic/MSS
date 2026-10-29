function [T,zeta] = DPperiods(vessel,display)
% DPperiods computes the DP periods and relative damping factors
%
%  [T,zeta] = DPperiods(vessel,display)
%
%  Inputs:
%     vessel    : MSS vessel structure
%     display   : 0 no display, 1 display results
%
%  Outputs:
%     T(1:6)    : Time constants in surge, sway, and yaw
%               : natural periods in heave, roll, and pitch
%     zeta(1:6) : relative damping factors in heave, roll, and pitch
%               : -1 for surge, sway and yaw
%
% Author:    Thor I. Fossen
% Date:      2005-09-26
% Revisions:
%   2026-03-30 : New formulas for periods/time constants based on vesselPeriods.m
%   2026-10-28 : Read viscous damping from vessel.powerBased.Bv.

if nargin == 1
    display = 0;
end

w   = vessel.freqs;
Nw = length(w);

MRB = vessel.MRB;
G   = reshape(vessel.C(:,:,Nw,1),6,6);
flagBV = isfield(vessel,'powerBased') && isfield(vessel.powerBased,'Bv');

for k = 1:Nw
    A(:,:,k) = reshape(vessel.A(:,:,k,1),6,6,1);

    if flagBV  % Viscous damping from the current vessel structure
        B(:,:,k) = reshape(vessel.B(:,:,k,1),6,6,1) + vessel.powerBased.Bv(:,:,1);
    else                     % No viscous damping
        B(:,:,k) = reshape(vessel.B(:,:,k,1),6,6,1);
    end

end

% *************************************************************************
% Compute periods/time constants
% *************************************************************************

% LF model DOFs 1,2,6 (uses first frequency as zero frequency)
M11 = MRB(1,1) + A(1,1,1);
M22 = MRB(2,2) + A(2,2,1);
M66 = MRB(6,6) + A(6,6,1);

if flagBV == 0  % if no visocus damping, use 20% of max Bii as LF estimate
    B11 = 0.2 * max(B(1,1,:));
    B22 = 0.2 * max(B(2,2,:));
    B66 = 0.2 * max(B(6,6,:));
else
    B11 = B(1,1,1);
    B22 = B(2,2,1);
    B66 = B(6,6,1);
end

% Heave, roll and pitch periods
[T345,zeta345] = vesselPeriods(w,vessel.MRB,A,B,G,'coupled',0);

% *************************************************************************
% Assign values to outputs
% *************************************************************************
T(1) = M11/B11;     % Time constants DOFs 1,2,6
T(2) = M22/B22;
T(6) = M66/B66;
T(3) = T345(1);     % Damped periods DOFs 3,4,5
T(4) = T345(2);
T(5) = T345(3);

% Damping factors
zeta(1) = -1;
zeta(2) = -1;
zeta(3) = zeta345(1);
zeta(4) = zeta345(2);
zeta(5) = zeta345(3);
zeta(6) = -1;

% *************************************************************************
% Display
% *************************************************************************
if display == 1
    disp('Time constant       Period    Damping')
    disp('------------------------------------------------')
    fprintf('Surge: %7.2f (s)\n',T(1))
    fprintf('Sway:  %7.2f (s)\n',T(2))
    fprintf('Heave: %18.2f (s) %6.3f\n',T(3),zeta(3))
    fprintf('Roll:  %18.2f (s) %6.3f\n',T(4),zeta(4))
    fprintf('Pitch: %18.2f (s) %6.3f\n',T(5),zeta(5))
    fprintf('Yaw:   %7.2f (s)\n',T(6))
    disp('------------------------------------------------')
end
