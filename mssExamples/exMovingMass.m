% exMovingMass Pitch control of an AUV using an internal moving mass.
% This example simulates the pitch and actuator loops for a 30.9 kg hull
% carrying a 1 kg mass that translates along the longitudinal axis. The
% vehicle is neutrally buoyant, with the nominal CG below the CB.
%
% The pitch controller sets a desired mass position. A faster PD servo
% drives the mass, with reference clipping and a cubic soft-wall force.
% The soft-wall force resists travel toward the end stops. It does not
% guarantee the physical travel limit for arbitrary commands or gains.
%
% The default pitch command changes from -30 deg to 10 deg at t = 60 s.
% The first command requires more mass travel than is available. These
% defaults reproduce the moving-mass example in Chapter 9.
%
% Positive pitch is nose-up. For an outer depth loop using NED depth,
% z_dot = -u*sin(theta) when heave velocity and roll are neglected.
% With u > 0, use the positive-error pitch command
%   theta_d = Kp_z*(z-z_d) + Ki_z*z_int;
%   z_int_dot = z-z_d;
% with positive gains. This example prescribes theta_d directly and
% simulates the inner pitch and mass-position loops.
%
% Usage:
%   exMovingMass
%
% Run in base MATLAB or GNU Octave. Parameters and simulation histories
% remain in the workspace. The example displays plots without exporting
% figure files.
%
% Reference:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK, Chapter 9.
%
% Author:    Thor I. Fossen
% Date:      2026-11-04
% Revisions:

clearvars; clc; close all

% Vehicle and moving-mass parameters
m_h = 30.9;              % Hull mass (kg)
m_p = 1.0;               % Moving mass (kg)
d_p = 0.1;               % Moving-mass damping (kg/s)
g = 9.81;                % Gravitational acceleration (m/s^2)
m = m_h + m_p;           % Total vehicle mass (kg)
W = m*g;                 % Weight (N)
M55 = 7.8273;            % Pitch inertia including added inertia (kg m^2)
zG = 0.02;               % Nominal CG below CB (m)

xmax = 0.30;             % Physical actuator travel limit (m)
xb = 0.22;               % Soft-wall boundary, xb < xmax (m)
kb = 5e3;                % Cubic soft-wall stiffness (N/m^3)

% Hydrostatic pitch stiffness and damping
Ktheta0 = zG*W;                  % Pitch stiffness (N m/rad)
wn0 = sqrt(Ktheta0/M55);         % Natural pitch frequency (rad/s)
zeta_rel = 0.15;                 % Relative pitch damping ratio
Dtheta = 2*zeta_rel*M55*wn0;      % Pitch damping (N m s/rad)

% Feedback-linearizing pitch controller
wn_th = 0.5;                    % Desired pitch bandwidth (rad/s)
zeta_th = 1.0;
Kp_th = M55*wn_th^2;
Kd_th = max(2*zeta_th*M55*wn_th-Dtheta,0);

% Moving-mass PD controller, faster than the pitch loop
wn_x = 1.0;                     % Mass-position bandwidth (rad/s)
zeta_x = 1.0;
Kp_x = m_p*wn_x^2;
Kd_x = max(2*zeta_x*m_p*wn_x-d_p,0);

% Simulation setup
h = 0.01;                       % Time step (s)
T = 120;                        % Duration (s)
t = 0:h:T;
N = numel(t);
lastSample = N;

theta = zeros(1,N);
q = zeros(1,N);
xP = zeros(1,N);
vP = zeros(1,N);
xPd = zeros(1,N);
tau_p = zeros(1,N);
tau_ca = zeros(1,N);

theta_d = deg2rad(-30)*ones(1,N);
theta_d(t >= T/2) = deg2rad(10);

for k = 1:N-1
    % Pitch-control effectiveness, with phi = 0
    b5 = -m_p*g*cos(theta(k));
    e_theta = atan2(sin(theta(k)-theta_d(k)), ...
                   cos(theta(k)-theta_d(k)));

    % Pitch controller: desired moving-mass position
    xPd(k) = (Ktheta0*sin(theta(k)) ...
              - Kp_th*e_theta - Kd_th*q(k))/b5;
    xPd(k) = max(min(xPd(k),xmax),-xmax);

    % Moving-mass servo and soft-wall force
    tau_ca(k) = -kb*sign(xP(k))*max(0,abs(xP(k))-xb)^3;
    tau_p(k) = -Kp_x*(xP(k)-xPd(k))-Kd_x*vP(k)+tau_ca(k);

    % Update velocities before positions
    vPdot = (tau_p(k)-d_p*vP(k))/m_p;
    vP(k+1) = vP(k)+h*vPdot;
    xP(k+1) = xP(k)+h*vP(k+1);

    theta_ddot = (b5*xP(k)-Dtheta*q(k) ...
                  - Ktheta0*sin(theta(k)))/M55;
    q(k+1) = q(k)+h*theta_ddot;
    theta(k+1) = theta(k)+h*q(k+1);

    % Stop if an edited case exceeds the physical travel limit
    if abs(xP(k+1)) > xmax
        warning('Moving mass exceeded physical limit at t = %.2f s',t(k+1));
        lastSample = k+1;
        break
    end
end

% Keep every history at the same length if the simulation stops early
t = t(1:lastSample);
theta = theta(1:lastSample);
theta_d = theta_d(1:lastSample);
q = q(1:lastSample);
xP = xP(1:lastSample);
vP = vP(1:lastSample);
xPd = xPd(1:lastSample);
tau_p = tau_p(1:lastSample);
tau_ca = tau_ca(1:lastSample);
N = lastSample;

xPd(end) = xPd(max(end-1,1));
tau_p(end) = tau_p(max(end-1,1));
tau_ca(end) = tau_ca(max(end-1,1));

% Plots
figure(1);
subplot(311)
plot(t,rad2deg(theta),'-.k','LineWidth',1.5); hold on
plot(t,rad2deg(theta_d),'-k','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('\theta (deg)')
xlim([0 t(end)])
ylim([-40 20])
legend('\theta','\theta_d')
title('Desired and actual pitch angles')
hold off

subplot(312)
plot(t,xP,'-.k','LineWidth',1.5); hold on
plot(t,xPd,'-k','LineWidth',1.5)
plot(t,xmax*ones(size(t)),'--k','LineWidth',1.0)
plot(t,-xmax*ones(size(t)),'--k','LineWidth',1.0)
grid on
xlabel('Time (s)')
ylabel('x_p (m)')
xlim([0 t(end)])
ylim([-0.4 0.4])
title('Moving-mass position')
legend('x_p','x_{p,d}','x_{max}','-x_{max}')
hold off

subplot(313)
plot(t,tau_p,'-k','LineWidth',1.5); hold on
plot(t,tau_ca,'-.k','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('Force (N)')
xlim([0 t(end)])
title('Moving-mass actuator force')
legend('\tau_p','\tau_{ca}')
hold off

set(findall(gcf,'type','axes'),'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',14);
set(findall(gcf,'type','legend'),'FontSize',12);
