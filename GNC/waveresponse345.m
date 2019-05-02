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
% T_0 = wave periode (s) corresponding to the wave frequency w_0 = 2 pi/T_0
% U = ship speed (m/s)
%
% Optionally ship data: ship = [zeta4, T4, GM_T, Cwp, Cb, delta]
% 
% zeta4      Relative damping factor in roll
% T4         Natural roll periode (s)
% GM_T       Transver metacentric height (m)
% Cb         Block coefficeint 
% Cw         Water plane area coefficient
%
% xdot = waveresponse345(t, x, a, beta, T_0, U, L, B, T, ship) allows the 
% user to specify the ship parameters. 
%
% xdot = waveresponse345(t, x, a, beta, T_0, U, L, B, T) uses the default 
% values  = [0.2, 6, 1, 0.65, 1.0].
%
% waveresponse345(-1, x, a, beta, T0, U, L, B, T,ship) solves the ODEs 
% and plots the heave, roll and pitch response for 20 seconds. 
% 
% Ex: waveresponse345(-1, zeros(6,1), 2, 230*pi/180, 10, 5, 82.8, 19.2, 6);
%
% Referances: 
% [1] J. Juncher Jensen, A. E. Mansour and A. S. Olsen. Estimation of 
% Ship Motions using Closed-form Expressions. Ocean Eng. 31, 2004, pp. 61-85
% 
% [2] O. M. Faltinsen. Sea Loads on Ships and Ocean Structures. Cambridge
% University Press, 1990.
% 
% Author:    Thor I. Fossen
% Date:      2018-07-21  Based on the method of Jensen et al. (2005)
% Revisions: 2019-05-xx  Added the method of Faltinen (1990)

% Default ship parameters
if nargin == 9
    ship = [0.2, 6, 1, 0.65, 1.0];
end

% Hydrodynamic parameters
zeta4   = ship(1);    % Relative damping factor in roll
T4      = ship(2);    % Natural roll periode (s)
GM_T    = ship(3);    % Transver metacentric height (m)
Cb      = ship(4);    % Block coefficeint 
Cw      = ship(5);    % Water plane coefficient

% Constants
g = 9.81;                 % acceleration of gravity (m/s^2)
rho = 1025;               % density of water (kg/m^3)

% Ship parameters
Aw = Cw * L * B;                % water plane area (m^2)
nabla = Cb * L * B * T;         % volume displacement (m^3) 
m = rho * nabla;                % mass
w_0 = 2 * pi / T_0;             % wave peak frequency (rad/s)
k = w_0^2/g;                    % wave number
w_e = w_0 - k * U * cos(beta);  % frequency of encounter
k_e = abs(k * cos(beta));       % effective wave number
sigma = k_e * L/2;
kappa = exp(-k_e * T);

% Heave and pitch models (Jensen et al., 2004)
alpha = w_e/w_0;
A = 2 * sin(k*B*alpha^2/2) * exp(-k*T*alpha^2); 
f = sqrt( (1-k*T)^2  + (A^2/(k*B*alpha^3))^2 );
F = kappa * f * sin(sigma)/sigma;
G = kappa * f * (1/sigma) * (6/L) * ( sin(sigma)/sigma - cos(sigma) );

% Heave and pitch models (Faltinsen, 1990)
A11 = 0.1 * m;              % 3-D added mass in surge
A11_2D = A11/T;             % 2-D added mass in surge
A33_2D = 0.8*rho*B*T;       % 2-D added mass in heave
A33 = A33_2D * L;           % 3-D added mass in heave         
M33 = m + A33;
R55 = 0.27 * L;
M55 = m * R55^2 + (1/12) * (A11_2D * T^3 + A33_2D * L^3);
C33 = rho * g * Aw;
F3 = g * L * (kappa*rho*a*B-1) * sin(sigma)/sigma;
F5 = (rho*g*a*2*B/k_e^2) * (kappa-0.8*k_e*T*exp(-k_e*T/2))*(sin(sigma)-sigma*cos(sigma));

% Natural frequency in Jensen et al. (2004) uses: w3 = w5 = sqrt(g/(2*T))
% Solution below uses spring stiffness and mass/added mass while relative
% damping ratio is based on Jensen et al. (2004)
w3 = sqrt(C33/M33);
w5 = w3;
zeta3 = (A^2/(B*alpha^3)) * sqrt(1/(8*k^3*T));
zeta5 = zeta3;

% Roll model (simplifed version of Jensen et al., 2004)
w4 = 2*pi/T4;                    % natural frequency
C44 = rho * g * nabla * GM_T;    % spring coeffient
M44 = C44/w4^2;                  % moment of inertia including added mass
B44 = 2 * zeta4 * w4 * M44;      % damping coefficient
M = sin(beta) * sqrt( B44 * rho*g^2/w_e );   % roll moment amplitude

% Output: xdot = f(x) 
if t == -1
    xdot = zeros(6,1);
else
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = w3^2 * (a*F*cos(w_e*t) - A^2/(k*B*alpha^3*w_0)*x(1) - x(4));
    xdot(5) = (1/M44) * (M*cos(w_e*t) - B44*x(2) - C44*x(5));
    xdot(6) = w5^2 * (a*G*sin(w_e*t) - A^2/(k*B*alpha^3*w_0)*x(3) - x(6));
end

% Optionally plotting. The solution of the ODE is valid for all frequencies 
% including the resonance where w_0 is equal to the natural frequency.
% URL: https://en.wikipedia.org/wiki/Harmonic_oscillator
if t == -1
      
    t = 0:0.1:20;    % time vector for plotting
        
    % Heave
    Z3 = sqrt( (2*w3*zeta3)^2 + (1/w_e^2)*(w3^2-w_e^2)^2 );
    eps3 = atan( 2*w_e*w3*zeta3/(w3^2-w_e^2) );
    z_response1 = (a*F*w3^2/(Z3*w_e)) * cos(w_e*t+eps3);
    z_response2 = (F3/(M33*Z3*w_e)) * cos(w_e*t+eps3);
    
    % Pitch
    Z5 = sqrt( (2*w5*zeta5)^2 + (1/w_e^2)*(w5^2-w_e^2)^2 );
    eps5 = atan( 2*w_e*w5*zeta5/(w5^2-w_e^2) );
    theta_response1 = (180/pi) * (a*G*w5^2/(Z5*w_e)) * sin(w_e*t+eps5);
    theta_response2 = (180/pi) * (F5/(M55*Z5*w_e))   * sin(w_e*t+eps5);    
    
    % Roll
    Z4 = sqrt( (2*w3*zeta4)^2 + (1/w_e^2)*(w4^2-w_e^2)^2 );
    eps4 = atan( 2*w_e*w4*zeta4/(w4^2-w_e^2) );
    roll_response = (180/pi) * ((M/C44)*w4^2/(Z4*w_e)) * cos(w_e*t+eps4);
    
    % Plots 
    figure(gcf)
    subplot(311)
    plot(t,z_response1,'-k',t,z_response2,'-.r','linewidth',2),xlabel('time (s)'), grid
    title(sprintf('Steady-state heave response (deg) for wave amplitude a = %2.1f m',a))
    legend('Jensen et al. (2004)','Faltinsen (1990)')
    subplot(312)
    plot(t,roll_response,'linewidth',2),xlabel('time (s)'), grid
    title(sprintf('Steady-state roll resopnse (deg) for wave amplitude a = %2.1f m',a))
    legend('Jensen et al. (2004)')
    subplot(313)
    plot(t,theta_response1,'-k',t,theta_response2,'-.r','linewidth',2),xlabel('time (s)'),grid
    title(sprintf('Steady-state pitch response (deg)for wave amplitude a = %2.1f m',a))
    legend('Jensen et al. (2004)','Faltinsen (1990)')
end


