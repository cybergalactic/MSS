function [T,zeta] = DPperiods(vessel,display)
%% DPperiods    
% Computes the following DP charachteristics:
%
%     >> [T,zeta] = DPperiods(vessel,display)
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
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%%
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
                 + reshape(vessel.B(:,:,k,1),6,6,1);
        flagBV = 1;
    else % no viscous damping
        B(:,:,k) = reshape(vessel.B(:,:,k,1),6,6,1);
        flagBV = 0;
    end
    
end

% *************************************************************************
% Compute periods/time constants for diagonal matrices 
% *************************************************************************

% estimate the natural frequencies in DOF 3,4,5
w_0 = sqrt(G(3,3)/(MRB(3,3)+A(3,3)));
w3 = natfrequency(vessel,3,w_0,1);
w_0 = sqrt(G(4,4)/(MRB(4,4)+A(4,4)));
w4 = natfrequency(vessel,4,w_0,1);
w_0 = sqrt(G(5,5)/(MRB(5,5)+A(5,5)));
w5 = natfrequency(vessel,5,w_0,1);

% LF model DOF 1,2,6
w_LF = 2*pi/150;
idx_126 = find(w>w_LF,1,'first');  % pick frequency >0 for period 150 s
  
M11 = MRB(1,1)+A(1,1,idx_126);
M22 = MRB(2,2)+A(2,2,idx_126);
M66 = MRB(6,6)+A(6,6,idx_126);

if flagBV == 0  % if no visocus damping, use 20% of max Bii as LF estimate
    B11 = 0.2*max(B(1,1,:));
    B22 = 0.2*max(B(2,2,:));
    B66 = 0.2*max(B(6,6,:));    
else
    B11 = B(1,1,idx_126);
    B22 = B(2,2,idx_126);
    B66 = B(6,6,idx_126);
end

Tdiag1 = [...
    M11/B11; 
    M22/B22;
    M66/B66 ];


[Tdiag1,idx1] = sort(Tdiag1);

% *************************************************************************
%% Compute eigenvalues, damping factors for coupled system
% *************************************************************************

% LF added mass
MA = A(:,:,idx_126);
MA(3,3) = interp1(w,reshape(A(3,3,:),1,Nw),w3);
MA(4,4) = interp1(w,reshape(A(4,4,:),1,Nw),w4);

% added mass at resonance freq
MA(5,5) = interp1(w,reshape(A(5,5,:),1,Nw),w5);

% total mass
M = MRB + MA;

% LF damping
N = B(:,:,idx_126);

if flagBV == 0  % if no visocus damping, use 20% of max Bii as LF estimate
    N(1,1) = 0.2*max(B(1,1,:));
    N(2,2) = 0.2*max(B(2,2,:));
    N(6,6) = 0.2*max(B(6,6,:));    
end
    
% damping at resonance freqs
N(3,3) = interp1(w,reshape(B(3,3,:),1,Nw),w3); 
N(4,4) = interp1(w,reshape(B(4,4,:),1,Nw),w4);
N(5,5) = interp1(w,reshape(B(5,5,:),1,Nw),w5);

if isfield(vessel,'roll') % only veres data
    N(4,4) = N(4,4) + interp1(w,reshape(vessel.roll.Bv44(:,length(vessel.velocities)),Nw,1),w4);
end

% system matrix
As = [...
    zeros(6,6) eye(6,6)
   -inv(M)*G -inv(M)*N ];

lam = eig(As);

k1 = 1;
k2 = 1;
Told = 0;
eps = 0.0001;

for i = 1:12

    wn  = abs(lam(i));
            
    if abs(imag(lam(i)) < eps)     % no spring (time constant)

        if real(lam(i)) < -eps  
            T  = 1/wn;         
            data126(k1,:) = T;
            k1 = k1 + 1;
        end        
    else                  % mass-damper-spring (period)
        T  = 2*pi/wn;          
        z = -real(lam(i))/wn;       
        
        if abs(T-Told) > eps
            data345(k2,:) = [T, z];
            k2 = k2+1;
            Told = T;
        end
    end

end

T126 = sort(data126);

[TT,idx] = sort(data345);
T345 = TT(:,1);
index = idx(:,1);

z345 = [...
    data345(index(1),2)
    data345(index(2),2)
    data345(index(3),2) ];

% *************************************************************************
%% Assign values to outputs
% *************************************************************************
T(1) = T126(idx1(idx1(1)));
T(2) = T126(idx1(idx1(2)));
T(6) = T126(idx1(idx1(3)));

% T(3) = T345(idx2(idx2(1)));
% T(4) = T345(idx2(idx2(2)));
% T(5) = T345(idx2(idx2(3)));
T(3) = 2*pi/w3;
T(4) = 2*pi/w4;
T(5) = 2*pi/w5;

zeta(1) = -1;
zeta(2) = -1;
zeta(6) = -1;

% coupled values in CG - must be updated to handle different resonance freqs.

% zeta(3) = z345(idx2(idx2(1)));
% zeta(4) = z345(idx2(idx2(2)));
% zeta(5) = z345(idx2(idx2(3)));

M33 = MRB(3,3) + interp1(w,reshape(A(3,3,:),1,Nw),w3);
M44 = MRB(4,4) + interp1(w,reshape(A(4,4,:),1,Nw),w4);
M55 = MRB(5,5) + interp1(w,reshape(A(5,5,:),1,Nw),w5);

N33 = interp1(w,reshape(B(3,3,:),1,Nw),w3);

if flagBV == 1  % for roll peack viscous value at w4 peak
    N44 = interp1(w,reshape(vessel.B(4,4,:,1),1,Nw),w4) + max(vessel.Bv(4,4,:));
else
    N44 = interp1(w,reshape(B(4,4,:),1,Nw),w4);
end

N55 = interp1(w,reshape(B(5,5,:),1,Nw),w5);

% 2*zeta*w = B/M  or  zeta = 0.5*(B/M)*(1/w)
zeta(3) = 0.5*(N33/M33)*(1/w3);
zeta(4) = 0.5*(N44/M44)*(1/w4);
zeta(5) = 0.5*(N55/M55)*(1/w5);

% *************************************************************************
%% Display
% *************************************************************************

if display == 1
    if flagBV == 1
        disp(' ')
        disp('ADDITIONAL SKIN FRICTION/ROLL DAMPING (U=0 m/s):')
        disp('Time constant       Period    Damping')
        disp('------------------------------------------------')
        disp(sprintf('Surge: %7.2f (s)',T(1)))
        disp(sprintf('Sway:  %7.2f (s)',T(2)))
        disp(sprintf('Heave: %18.2f (s) %6.3f',T(3),zeta(3)))
        disp(sprintf('Roll:  %18.2f (s) %6.3f',T(4),zeta(4)))
        disp(sprintf('Pitch: %18.2f (s) %6.3f',T(5),zeta(5)))
        disp(sprintf('Yaw:   %7.2f (s)',T(6)))
        disp('------------------------------------------------')
    else
        disp(' ')
        disp('DAMPING FROM HYDRODYNAMIC CODE (U=0 m/s):')
        disp('Natural frequency    Periods   Damping')
        disp('-----------------------------------------')
        disp(sprintf('Heave: %3.2f (rad/s) %6.2f (s) %6.3f',2*pi/T(3),T(3),zeta(3)))
        disp(sprintf('Roll:  %3.2f (rad/s) %6.2f (s) %6.3f',2*pi/T(4),T(4),zeta(4)))
        disp(sprintf('Pitch: %3.2f (rad/s) %6.2f (s) %6.3f',2*pi/T(5),T(5),zeta(5)))
        disp('-----------------------------------------')
    end
end



