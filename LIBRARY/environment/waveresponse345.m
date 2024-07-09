function waveresponse345(a, beta,T_0, zeta4,T4,GMT, Cb, U, L, B, T)
% waveresponse345(a, beta,T_0, zeta4,T4,GMT, Cb, U, L, B, T) computes and 
% plots the steady-state heave, roll and pitch responses for a ship in 
% regular waves using closed-form formulae. 
%
% Ex: waveresponse345(2, 75*pi/180, 10, 0.2, 6, 1, 0.65, 5, 82.8, 19.2, 6)
%
% Inputs:
%
% a      Wave amplitude (m)
% beta   Wave direction (rad) where beta = pi (i.e. 180 deg) is head seas
% T_0    Wave periode (s) corresponding to the wave frequency w_0 = 2 pi/T_0
% zeta4  Relative damping factor in roll
% T4     Natural roll periode (s)
% GMT    Transver metacentric height (m)
% Cb     Block coefficeint 
% U      Ship speed (m/s)
% L      Length (m)
% B      Breadth (m)
% T      Draught (m)
% 
%
% Reference: 
%   J. Juncher Jensen, A. E. Mansour and A. S. Olsen. Estimation of 
%   Ship Motions using Closed-form Expressions. Ocean Eng. 31, 2004, pp.61-85
% 
% Author:    Thor I. Fossen
% Date:      2018-07-21  
% Revisions: 
%   2019-05-04  Bug fixes
%   2019-05-11  Updated formula for wn, removed bugs

% Constants
g = 9.81;                 % acceleration of gravity (m/s^2)
rho = 1025;               % density of water (kg/m^3)

% Ship parameters
nabla = Cb * L * B * T;        % volume displacement (m^3) 
w_0 = 2 * pi / T_0;            % wave peak frequency (rad/s)
k = w_0^2/g;                   % wave number
w_e = w_0 - k * U * cos(beta); % frequency of encounter
k_e = abs(k * cos(beta));      % effective wave number
sigma = k_e * L/2;
kappa = exp(-k_e * T);

% Heave and pitch models (Jensen et al., 2004)
alpha = w_e/w_0;
A = 2 * sin(k*B*alpha^2/2) * exp(-k*T*alpha^2); 
f = sqrt( (1-k*T)^2  + (A^2/(k*B*alpha^3))^2 );
F = kappa * f * sin(sigma)/sigma;
G = kappa * f * (6/L) * (1/sigma) * ( sin(sigma)/sigma - cos(sigma) );

% Natural frequency in Jensen et al. (2004) uses: w3 = w5 = sqrt(g/(2*T))
% Solution below uses spring stiffness and mass/added mass while relative
% damping ratio is based on Jensen et al. (2004)
wn  = sqrt(g/(2*T));
zeta = (A^2/(B*alpha^3)) * sqrt(1/(8*k^3*T));

% Roll model (simplifed version of Jensen et al., 2004)
w4 = 2*pi/T4;                    % natural frequency
C44 = rho * g * nabla * GMT;     % spring coeffient
M44 = C44/w4^2;                  % moment of inertia including added mass
B44 = 2 * zeta4 * w4 * M44;      % damping coefficient
M = sin(beta) * sqrt( B44 * rho*g^2/w_e );   % roll moment amplitude

% The solution of the ODE is valid for all frequencies 
% including the resonance where w_0 is equal to the natural frequency.
% URL: https://en.wikipedia.org/wiki/Harmonic_oscillator
t = 0:0.1:20;    % time vector for plotting
        
% Heave
Z3 = sqrt( (2*wn*zeta)^2 + (1/w_e^2)*(wn^2-w_e^2)^2 );
eps3 = atan( 2*w_e*wn*zeta/(wn^2-w_e^2) );
z = (a*F*wn^2/(Z3*w_e)) * cos(w_e*t+eps3);
    
% Pitch
Z5 = sqrt( (2*wn*zeta)^2 + (1/w_e^2)*(wn^2-w_e^2)^2 );
eps5 = atan( 2*w_e*wn*zeta/(wn^2-w_e^2) );
theta = (180/pi) * (a*G*wn^2/(Z5*w_e)) * sin(w_e*t+eps5);
    
% Roll
Z4 = sqrt( (2*w4*zeta4)^2 + (1/w_e^2)*(w4^2-w_e^2)^2 );
eps4 = atan( 2*w_e*w4*zeta4/(w4^2-w_e^2) );
phi = (180/pi) * ((M/C44)*w4^2/(Z4*w_e)) * cos(w_e*t+eps4);
    
%% Plots 
clf
figure(gcf)
hold on
plot(t,z,'-k','linewidth',2)
plot(t,phi,':r','linewidth',2)
plot(t,theta,'-.b','linewidth',2)
hold off
title(sprintf('Steady-state responses for a = %2.1f m and beta = %2.1f deg',a,(180/pi)*beta))
legend('Heave (m)','Roll (deg)','Pitch (deg)')
xlabel('time (s)')
grid



