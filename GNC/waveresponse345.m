function xdot = waveresponse345(t, x, a, beta,T_0, U, L, B, T, ship)
% WAVERESPONSE345 computes the wave-induced ship motions for a regular wave
% in heave, roll and pitch using closed-form formulae. The function returns 
% the time derivative xdot of the state vector x = [z, phi, theta, w, p, r]'.
%
% Inputs:
%
% t = time (s)
% x = [z (m), phi (rad), theta (rad), w (m/s), p (rad/s), r (rad(s)]'
% a = wave amplitude (m)
% beta = wave direction (rad) where beta = pi (i.e. 180 deg) is head seas
% T0 = wave periode (s) corresponding to the wave frequency w0 = 2 pi/T0
% U = ship speed (m/s)
%
% Optionally ship data: ship = [zeta_roll, T_roll, GM_T, Cwp, Cb, delta]
% 
% zeta_roll  Relative damping factor in roll
% T_roll     Natural roll periode (s)
% GM_T       Transver metacentric height (m)
% Cb         Block coefficeint 
% delta      0 < delta = bow/L < Cwp, see Jensen et al. 2014 for details
%
% xdot = waveresponse345(t, x, a, beta, T0, U, L, B, T, ship) allows the 
% user to specify the ship parameters. 
%
% xdot = waveresponse345(t, x, a, beta, T0, U, L, B, T) uses the default 
% values  = [0.2, 6, 1, 0.65, 1/3].
%
% waveresponse345(-1, x, a, beta, T0, U, L, B, T,ship) solves the ODEs 
% and plots the heave, roll and pitch response for 20 seconds. 
% 
% Ex: waveresponse345(-1, zeros(6,1), 2, 230*pi/180, 10, 5, 40, 6, 3);
%
% Ref: J. Juncher Jensen, A. E. Mansour and A. S. Olsen. Estimation of 
% ship motions using closed-form expressions
% Ocean Engineering 31, 2004, pp. 61-85
%
% Author:    Thor I. Fossen
% Date:      21 July 2018
% Revisions: 

if nargin == 9
    ship = [0.2, 6, 1, 0.65, 1/3];
end

% Constants
g = 9.81;                 % acceleration of gravity (m/s^2)
rho = 1025;               % density of water (kg/m^3)

% Hydrodynamic parameters
zeta_roll   = ship(1);    % Relative damping factor in roll
T_roll      = ship(2);    % Natural roll periode (s)
GM_T        = ship(3);    % Transver metacentric height (m)
Cb          = ship(4);    % Block coefficeint 
delta       = ship(5);    % 0 < delta = bow/L < Cwp

if ((delta > L*B) | (delta <= 0))
    error('0 < delta = bow/L < L*B');
end

nabla = Cb * L * B * T;         % volume displacement (m^3) 
w_0 = 2 * pi / T_0;             % Wave peak frequency (rad/s)
k = w_0^2/g;                    % wave number
w_e = w_0 - k * U * cos(beta);  % frequency of encounter

% Heave and pitch models
alpha = w_e/w_0;
k_e = abs(k * cos(beta));
KL = k_e * L/2;

kappa = exp(-k_e * T);
A = 2 * sin(k*B*alpha^2/2) * exp(-k*T*alpha^2); 
f = sqrt( (1-k*T)^2  + (A^2/(k*B*alpha^3))^2 );
F = kappa * f * (1/KL) * sin(KL);
G = kappa * f * (1/KL)^2 * (6/L) * ( sin(KL) - KL*cos(KL) );

% Roll model
w_roll = 2*pi / T_roll;              % natural frequency
C44 = rho * g * nabla * GM_T;        % spring coeffient
M44 = C44/w_roll^2;                  % moment of inertia including added mass
B44 = 2 * zeta_roll * w_roll * M44;  % damping coefficient
M = sin(beta) * sqrt( B44 * rho*g^2/w_e );   % roll moment amplitude

% Outputs
if t == -1
    xdot = zeros(6,1);
else
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = g/(2*T) * (a*F*cos(w_e*t) - A^2/(k*B*alpha^3*w_0)*x(1) - x(4));
    xdot(5) = 1/M44   * (M*cos(w_e*t) - B44*x(2) - C44*x(5));
    xdot(6) = g/(2*T) * (a*G*sin(w_e*t) - A^2/(k*B*alpha^3*w_0)*x(3) - x(6));
end

% Optionally plotting. The solution of the ODE is valid for all frequencies 
% including the resonance where w_0 is equal to the natural frequency.
% URL: https://en.wikipedia.org/wiki/Harmonic_oscillator

if t == -1
      
    t = 0:0.1:20;    % time vector for plotting
        
    % Heave: F0  = a F m w3^2
    w3 = sqrt( g/(2*T) );
    zeta3 = (1/(2*w3)) * (A^2/(k*B*alpha^3*w_0))*(g/(2*T));
    Z3 = sqrt( (2*w3*zeta3)^2 + (1/w_e^2)*(w3^2-w_e^2)^2 );
    eps3 = atan( 2*w_e*w3*zeta3/(w3^2-w_e^2) );
    z_response = (a*F*w3^2/(Z3*w_e)) * cos(w_e*t+eps3);
    
    % Pitch: F0  = a G m w5^2
    w5 = sqrt( g/(2*T) );
    zeta5 = (1/(2*w5)) * (A^2/(k*B*alpha^3*w_0))*(g/(2*T));
    Z5 = sqrt( (2*w5*zeta5)^2 + (1/w_e^2)*(w5^2-w_e^2)^2 );
    eps5 = atan( 2*w_e*w5*zeta5/(w5^2-w_e^2) );
    theta_response = (180/pi) * (a*G*w5^2/(Z5*w_e)) * sin(w_e*t+eps5);
    
    % Roll  F0 = (M/C44) m w4^2
    w4 = 2*pi/T_roll;
    zeta4 = zeta_roll;
    Z4 = sqrt( (2*w3*zeta4)^2 + (1/w_e^2)*(w4^2-w_e^2)^2 );
    eps4 = atan( 2*w_e*w4*zeta4/(w4^2-w_e^2) );
    roll_response = (180/pi) * ((M/C44)*w4^2/(Z4*w_e)) * cos(w_e*t+eps4);
    
    % Plots 
    figure(gcf)
    subplot(311)
    plot(t,z_response,'linewidth',2),title('Heave (m)'),xlabel('time (s)'), grid
    subplot(312)
    plot(t,roll_response,'linewidth',2),title('Roll (deg)'),xlabel('time (s)'), grid
    subplot(313)
    plot(t,theta_response,'linewidth',2),title('Pitch (deg)'),xlabel('time (s)'),grid
    
end


