% exOtterEKF Five-state EKF and PID course control for the Otter USV.
% This script simulates the Maritime Robotics Otter Uncrewed Surface 
% Vehicle (USV), which is controlled by two propellers. SOG, COG, and 
% course rate are estimated using an EKF, and the course autopilot is 
% utilizing a PID pole-placement algorithm.
% This example uses the MSS-Capytaine otter.m model.
% The thrust schedule gives transit near 0.5 m/s and a low-speed phase
% near 0.2 m/s to illustrate the limits of course estimation and control.
% Use SIMotter for path following and optional wave excitation.
%
% NOTE that the estimates of SOG and COG start to osciallate when the USV
% speed is smaller than or close to 0.2 m/s. In this case, the course angle is
% not well-defined nor observable by the EKF.
%
% Author:    Thor I. Fossen
% Date:      2019-07-18
% Revisions: 
%   2026-11-01 - Use the MSS-Capytaine Otter model. 

clearvars;
clear EKF_5states   % Clear persistent states in EKF_5states.m

%% Simulation data
fHz = 50;           % Sampling frequency [Hz]
h  = 1/fHz;         % Sampling time [s]
Z = 5;              % GNSS updates every 5 samples (10 Hz)
T_final = 100;	    % Final simulation time [s]

% Initial values for x = [ u v w p q r x y z phi theta psi ]'
x = zeros(12,1); x(1) = 0;	   

% EKF covariance matrices and initial states
Qd = 1e6 * diag([ 1 1]);
Rd = 1 * diag([ 1 1 ]);
x_prd = zeros(5,1);

% Otter USV input matrix
[~,~,M, B_prop] = otter();
% The /2 converts the combined twin-propeller coefficients to per-propeller
% values. These coefficients must match CRAFT/USV/models/otter.m.
k_pos = 0.02216/2;
k_neg = 0.01289/2;
% B_prop maps n.*abs(n) through k_pos. Dividing by k_pos gives the
% geometry-only matrix mapping individual propeller thrusts to [tau_X; tau_N].
B_thrust = B_prop / k_pos;
B_thrust_inv = invQR(B_thrust);

% PID heading autopilot parameters (Nomoto model: M(6,6) = T/K)
T = 1;                           % Nomoto time constant
K = T / M(6,6);                  % Nomoto gain constant

wn = 0.85;                        % Closed-loop natural frequency (rad/s)
zeta = 1.2;                       % Closed-loop relative damping factor (-)

Kp = M(6,6) * wn^2;                     % Proportional gain
Kd = M(6,6) * (2 * zeta * wn - 1/T);    % Derivative gain
Td = Kd / Kp;                           % Derivative time constant
Ti = 10 / wn;                           % Integral time constant
e_chi = 0;                              % Initial integral state

% Reference model parameters
wn_d = 0.3;                      % Natural frequency (rad/s)
zeta_d = 1.0;                    % Relative damping factor (-)
r_max = deg2rad(10.0);           % Maximum turning rate (rad/s)

chi_d = 0;
omega_d = 0;
a_d = 0;

% Propeller dynamics
Tn = 0.1;                        % Propeller time constant      
n = [0 0]';                      % n = [ n_left n_right ]'

% Load condition
mp = 25;                         % Payload mass (kg), maximum value 45 kg
rp = [0.05 0 -0.35]';            % Location of payload (m)

% Ocean current
V_c = 0;               % Current speed (m/s)
beta_c = deg2rad(30);  % Current direction (rad)

%% MAIN LOOP
t = 0:h:T_final;                 % Time vector
simdata = zeros(length(t),19);   % Table for simulation data
for i=1:length(t)

   % GNSS positions at sample k. x_prd holds the prior estimate at k.
   x_N = x(7); y_E = x(8);

   % Setpoints
   if t(i) < 40
       chi_ref = deg2rad(20);
       tau_X = 8.5;       % Transit speed approximately 0.5 m/s
   elseif t(i) > 20 && t(i) < 30
       chi_ref = deg2rad(0); 
       tau_X = 8.5;       % Transit speed approximately 0.5 m/s
   else
       chi_ref = deg2rad(0); 
       tau_X = 2.8;       % Low-speed phase to illustrate COG limitations
   end
   
   % Course autopilot (PID pole placement with reference model)
   course_error = ssa(x_prd(4) - chi_d);
   tau_N = (T/K) * a_d + (1/K) * omega_d -... 
        Kp * ( course_error +...
        Td * (x_prd(5) - omega_d) + (1/Ti) * e_chi );
   thrust_c = B_thrust_inv * [tau_X; tau_N];
   n_c = zeros(2,1);
   for j = 1:2
       if thrust_c(j) >= 0
           n_c(j) = sqrt(thrust_c(j) / k_pos);
       else
           n_c(j) = -sqrt(-thrust_c(j) / k_neg);
       end
   end
   
   % n_c = [80 60]';          % n = [ n_left n_right ]' 
   
   % Reference model propagation
   [chi_d, omega_d, a_d] = refModel(chi_d,omega_d,a_d,chi_ref, ...
      r_max, zeta_d,wn_d,h,1);
  
   % Store simulation data in a table   
   simdata(i,:) = [x' x_prd' n'];    
   
   % RK4 method x(k+1)
   x = rk4(@otter,h,x,n,mp,rp,V_c,beta_c);
   
   % Euler's method (k+1)
   e_chi = e_chi + h * course_error;
   n = n - h / Tn * (n - n_c);  

   % Correct with GNSS[k] and return the prior estimate for sample k+1.
   x_prd = EKF_5states(x_N,y_E,h,Z,'NED',Qd,Rd);

end

%% PLOTS
nu   = simdata(:,1:6); 
eta  = simdata(:,7:12); 

psi = eta(:,6);                  % Heading angle
velocity_n = zeros(length(t),3);
for i = 1:length(t)
    Rnb = Rzyx(eta(i,4),eta(i,5),eta(i,6));
    velocity_n(i,:) = (Rnb * nu(i,1:3)')';
end
U_g = hypot(velocity_n(:,1),velocity_n(:,2));
chi = atan2(velocity_n(:,2),velocity_n(:,1));
chi(U_g == 0) = NaN;             % COG is undefined at zero horizontal speed

x_hat = simdata(:,13);
y_hat = simdata(:,14);
U_g_hat = simdata(:,15);
chi_hat = simdata(:,16);
omega_hat = simdata(:,17);

n1 = simdata(:,18);
n2 = simdata(:,19);

clf
figure(1)
subplot(111),plot(t,n1,'g',t,n2,'k'),grid
xlabel('time (s)')
legend('n_1 left propeller','n_2 right propeller')
title('Propeller revolutions (rad/s)')
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',11)
set(findall(gcf,'type','line'),'linewidth',1.5)

figure(2)
subplot(621),plot(t,nu(:,1))
xlabel('time (s)'),title('Surge velocity (m/s)'),legend('True'),grid
subplot(623),plot(t,nu(:,2))
xlabel('time (s)'),title('Sway velocity (m/s)'),legend('True'),grid
subplot(625),plot(t,nu(:,3))
xlabel('time (s)'),title('Heave velocity (m/s)'),legend('True'),grid
subplot(627),plot(t,rad2deg(nu(:,4)))
xlabel('time (s)'),title('Roll rate (deg/s)'),legend('True'),grid
subplot(629),plot(t,rad2deg(nu(:,5)))
xlabel('time (s)'),title('Pitch rate (deg/s)'),legend('True'),grid
subplot(6,2,11),plot(t,rad2deg(nu(:,6)))
xlabel('time (s)'),title('Yaw rate (deg/s)'),legend('True'),grid

subplot(622)
plot(t,eta(:,1),t,x_hat); 
title('North position (m)'),legend('True','Estimate'), grid
subplot(624)
plot(t,eta(:,2),t,y_hat); 
title('East position (m)'),legend('True','Estimate'), grid
subplot(626),plot(t,eta(:,3))
xlabel('time (s)'),title('Down position (m)'),legend('True'),grid
subplot(628),plot(t,rad2deg(eta(:,4)))
xlabel('time (s)'),title('Roll angle (deg)'),legend('True'),grid
subplot(6,2,10),plot(t,rad2deg(eta(:,5)))
xlabel('time (s)'),title('Pitch angle (deg)'),legend('True'),grid
subplot(6,2,12),plot(t,rad2deg(psi))
xlabel('time (s)'),title('Yaw angle (deg)')

set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',11)
set(findall(gcf,'type','line'),'linewidth',1.5)

figure(3)
set(gcf,'Units','inches','Position',[1 1 6.5 4.4], ...
    'PaperUnits','inches','PaperPosition',[0 0 6.5 4.4], ...
    'PaperSize',[6.5 4.4],'PaperPositionMode','manual');
subplot(221),plot(eta(:,2),eta(:,1),'-.',y_hat,x_hat,'-'); 
title('North-East positions (m)'),legend('True','Estimate','Location','northwest'),axis 'equal'; grid
xlabel('East position (m)'),ylabel('North position (m)')
north_min = min([eta(:,1); x_hat]);
north_max = max([eta(:,1); x_hat]);
north_margin = 0.05 * max(1,north_max-north_min);
ylim([north_min-north_margin north_max+north_margin]);
subplot(222),plot(t,U_g,'-.',t,U_g_hat,'-');
xlabel('time (s)'),title('SOG (m/s)'),legend('True','Estimate'),grid
subplot(223),plot(t,rad2deg(omega_hat),'r')
xlabel('time (s)'),title('Course rate (deg/s)')
legend('Estimate','Location','northwest'),grid
subplot(224),plot(t,rad2deg(chi),'-.',t,rad2deg(chi_hat),'-')
xlabel('time (s)'),title('COG (deg)')
% Show the commanded course range; the initial direction near rest may lie outside it.
ylim([-20 25]);
legend('True','Estimate','Location','southwest'),grid

set(findall(gcf,'type','text'),'FontSize',10)
set(findall(gcf,'type','legend'),'FontSize',9)
set(findall(gcf,'type','line'),'linewidth',1.0)

% Font sizes for the four-panel EKF comparison.
set(findall(figure(3),'Type','axes'),'FontSize',9);
axes_handles = findall(figure(3),'Type','axes');
for i = 1:length(axes_handles)
    set(get(axes_handles(i),'Title'),'FontSize',11,'FontWeight','normal');
end
