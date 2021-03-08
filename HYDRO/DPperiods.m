function [T,zeta] = DPperiods(vessel,display)
% DPperiods computes the DP periods and relative damping factors
%
%  [T,zeta] = DPperiods(vessel,display)
%
%  Inputs:
%     vessel    : MSS vessel structure
%     display   : (optionally) 0 no display, 1 display results to screen
%
%  Outputs:
%     T(1:6)    : Time constants in surge, sway, and yaw
%               : natural periods in heave, roll, and pitch 
%     zeta(1:6) : relative damping factors in heave, roll, and pitch 
%               : -1 for surge, sway and yaw
%
% Author:    Thor I. Fossen
% Date:      2005-09-26
% Revisions: 2008-05-09  First version
%            2009-09-11  Using constant viscous damping Bv
%            2013-07-06  Using new initial values w_0
%            2019-04-25  Fixed bug for test on imag and real eigenvalues
% Revisions: 2021-03-07  Major improvements, new formulas for periods/time constants

if nargin == 1
    display = 0;
end

w   = vessel.freqs;
Nw = length(w);
 
MRB = vessel.MRB;
G   = reshape(vessel.C(:,:,Nw,1),6,6);

for k = 1:Nw
    A(:,:,k) = reshape(vessel.A(:,:,k,1),6,6,1);
       
    if isfield(vessel,'Bv')  % viscous damping
        B(:,:,k) = reshape(vessel.B(:,:,k,1),6,6,1)... 
                 + vessel.Bv(:,:,1);
        flagBV = 1;
    else                     % no viscous damping
        B(:,:,k) = reshape(vessel.B(:,:,k,1),6,6,1);
        flagBV = 0;
    end
    
end

% *************************************************************************
%% Compute periods/time constants for diagonal matrices 
% *************************************************************************

% estimate the natural frequencies in DOFs 3,4,5
w_0 = sqrt(G(3,3)/(MRB(3,3)+A(3,3))); w3 = natfrequency(vessel,3,w_0,1);
w_0 = sqrt(G(4,4)/(MRB(4,4)+A(4,4))); w4 = natfrequency(vessel,4,w_0,1);
w_0 = sqrt(G(5,5)/(MRB(5,5)+A(5,5))); w5 = natfrequency(vessel,5,w_0,1);

% LF model DOFs 1,2,6 (uses first frequency as zero frequency)
M11 = MRB(1,1) + A(1,1,1);
M22 = MRB(2,2) + A(2,2,1);
M66 = MRB(6,6) + A(6,6,1);

if flagBV == 0  % if no visocus damping, use 20% of max Bii as LF estimate
    B11 = 0.2 * max(B(1,1,:));
    B22 = 0.2 * max(B(2,2,:));
    B66 = 0.2 * max(B(6,6,:));    
else
    B11 = B(1,1,1) ;
    B22 = B(2,2,1);
    B66 = B(6,6,1);
end

% *************************************************************************
%% Compute eigenvalues, damping factors for coupled system
% *************************************************************************

% LF added mass
MA = A(:,:,1);
MA(3,3) = interp1(w,reshape(A(3,3,:),1,Nw),w3);
MA(4,4) = interp1(w,reshape(A(4,4,:),1,Nw),w4);
MA(5,5) = interp1(w,reshape(A(5,5,:),1,Nw),w5);

% total mass
M = MRB + MA;

% Damping at resonant frequencies
N = B(:,:,1);
N(3,3) = interp1(w,reshape(B(3,3,:),1,Nw),w3); 
N(4,4) = interp1(w,reshape(B(4,4,:),1,Nw),w4);
N(5,5) = interp1(w,reshape(B(5,5,:),1,Nw),w5);

if isfield(vessel,'roll') % only veres data
    N(4,4) = N(4,4) + interp1(w,reshape(vessel.roll.Bv44(:,length(vessel.velocities)),Nw,1),w4);
end

% *************************************************************************
%% Assign values to outputs
% *************************************************************************
T(1) = M11/B11;     % Time constants DOFs 1,2,6
T(2) = M22/B22;
T(6) = M66/B66;

T(3) = 2*pi/w3;     % Natural periods DOFs 3,4,5
T(4) = 2*pi/w4;
T(5) = 2*pi/w5;

% coupled values in CG - must be updated to handle different resonance freqs.
M33 = MRB(3,3) + interp1(w,reshape(A(3,3,:),1,Nw),w3);
M44 = MRB(4,4) + interp1(w,reshape(A(4,4,:),1,Nw),w4);
M55 = MRB(5,5) + interp1(w,reshape(A(5,5,:),1,Nw),w5);

N33 = interp1(w,reshape(B(3,3,:),1,Nw),w3);

if flagBV == 1  % for roll peack viscous value at w4 
    N44 = interp1(w,reshape(vessel.B(4,4,:,1),1,Nw),w4) + max(vessel.Bv(4,4,:));
else
    N44 = interp1(w,reshape(B(4,4,:),1,Nw),w4);
end

N55 = interp1(w,reshape(B(5,5,:),1,Nw),w5);

% Damping factors
zeta(1) = -1;
zeta(2) = -1;
zeta(3) = 0.5 * (N33/M33) * (1/w3);
zeta(4) = 0.5 * (N44/M44) * (1/w4);
zeta(5) = 0.5 * (N55/M55) * (1/w5);
zeta(6) = -1;

% *************************************************************************
%% Display
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



