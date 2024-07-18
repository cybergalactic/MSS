function SIMrig()
% SIMrig is compatibel with MATLAB and GNU Octave (www.octave.org). This
% script simulates a semisubmersible with two pontoons and four columns. 
% The model matrices MRB, MA, D and G are generated using the utility script 
% 'rig.m', which calls the data file 'data_rig.m'. 
%
% Main characteristics:
%   m = 27 162 500 kg (mass)
%   H_p = 13.0 m (pontoon height)
%   B_p = 12.0 m (pontoon width)
%   L_p = 84.6 m (pontoon length)
%   r_CG = [-0.5 0 -19.5]', CG with respect to the CO
%
% Author:     Thor I. Fossen
% Date:       2019-03-12
% Revisions:  
%   2020-10-22 : Tuning of the PID controller.
%   2020-08-04 : Major updates.
%   2024-03-27 : DP control law replaced by PIDnonlinearMIMO.m. 
%   2024-04-19 : Enhanced compatibility with GNU Octave.

clearvars;
close all;
clear PIDnonlinearMIMO;        % Clear persistent integrator state

T_final = 200;	               % Final simulation time (s)
h = 0.1;                       % Sampling time (s)

% Environmnetal data
nu_c = [-1.0 -0.2 0 0 0 0]';    % ocean current velocities

% Nonlinear PID controller  
psi_d = deg2rad(10);            % Heading angle step
eta_ref = [0 0 psi_d]';         % DP system setpoints
wn = 0.2;                       % Closed-loop natural frequency (rad/s)
zeta = 1.0;                     % Closed-loop relative damping factor  
T_f = 5;                        % Low-pass filter time constant (s)

% Initial values
eta = [0 0 0 deg2rad(10) deg2rad(5) 0]';   % Generalized position vector
nu  = zeros(6,1);                          % Generalized velocity vector

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% Semisub model matrices obtained by running 'rig.m'
MRB = 1.0e+10 * [
    0.0027         0         0         0   -0.0530         0
    0    0.0027         0    0.0530         0   -0.0014
    0         0    0.0027         0    0.0014         0
    0    0.0530         0    3.4775         0   -0.0265
    -0.0530         0    0.0014         0    3.8150         0
    0   -0.0014         0   -0.0265         0    3.7192];

MA = 1.0e+10 * [
    0.0017         0         0         0   -0.0255         0
    0    0.0042         0    0.0365         0         0
    0         0    0.0021         0         0         0
    0    0.0365         0    1.3416         0         0
    -0.0255         0         0         0    2.2267         0
    0         0         0         0         0    3.2049];

D = 1.0e+09 * [
    0.0004         0         0         0   -0.0085         0
    0    0.0003         0    0.0067         0   -0.0002
    0         0    0.0034         0    0.0017         0
    0    0.0067         0    4.8841         0   -0.0034
    -0.0085         0    0.0017         0    7.1383         0
    0   -0.0002         0   -0.0034         0    0.8656];

G = 1.0e+10 * [
    0         0         0         0         0         0
    0         0         0         0         0         0
    0         0    0.0006         0         0         0
    0         0         0    1.4296         0         0
    0         0         0         0    2.6212         0
    0         0         0         0         0         0 ];

M = MRB + MA;    
Minv = inv(M);  

% Display simulation options
displayControlMethod();

%% MAIN LOOP 
simdata = zeros(nTimeSteps, 12); % Pre-allocate matrix for efficiency

for i=1:nTimeSteps
    
    % Measurements
    eta(1) = eta(1) + 0.01 * randn;
    eta(2) = eta(2) + 0.01 * randn;
    eta(6) = eta(6) + 0.0001 * randn;

    % DP control law
    if t(i) > 100
        eta_ref = [0 0 -psi_d]';
    end

    tau_thr = PIDnonlinearMIMO(eta,nu,eta_ref,M,wn,zeta,T_f,h);

    % Equations of motion
    nudot = Minv * ( tau_thr - D * (nu-nu_c) - G * eta );

    % Kinematics
    J = eulerang(eta(4),eta(5),eta(6));

    % Store data for presentation
    simdata(i,:) = [eta', nu'];

    % Euler's integration methods (k+1), (Fossen 2021, Eq. B27-B28)
    nu = nu + h * nudot;          % Forward Euler
    eta = eta + h * J * nu;       % Backward Euler

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]
legendLocation = 'best';
if isoctave; legendLocation = 'northeast'; end

x     = simdata(:,1); 
y     = simdata(:,2);          
z     = simdata(:,3);   
phi   = rad2deg(simdata(:,4)); 
theta = rad2deg(simdata(:,5)); 
psi   = rad2deg(simdata(:,6));
u     = simdata(:,7);
v     = simdata(:,8);
w     = simdata(:,9);
p     = rad2deg(simdata(:,10));
q     = rad2deg(simdata(:,11));
r     = rad2deg(simdata(:,12)); 

figure(1); 
if ~isoctave; set(gcf,'Position',[1, 1, scrSz(3)/3, scrSz(4)]); end
subplot(211)
plot(y,x)
axis('equal')
grid
xlabel('East')
ylabel('North')
title('Semisub xy-plot (m)')
subplot(212)
plot(t,x,t,y),title('Semisub positions (m)'),
legend('x position','y position','Location',legendLocation)
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend','Location',legendLocation),'FontSize',14)

figure(2); figure(gcf)
subplot(311)
plot(t,phi),xlabel('time (s)')
title('Roll angle \phi (deg)')
grid
subplot(312)
plot(t,theta),xlabel('time (s)')
title('Pitch angle \theta (deg)')
grid
subplot(313)
plot(t,psi)
hold on
plot( [0, t(end)],-rad2deg([psi_d, psi_d]),'c')
plot( [0, t(end)], rad2deg([psi_d, psi_d]),'c')
hold off
xlabel('time (s)')
title('Yaw angle \psi (deg)')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend','Location',legendLocation),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Pontion length', '84.6 m', ...
    'Pontoon width', '12.0 m', ...
    'Pontoon height', '13.0 m', ...
    'Width of box-shaped leg', '16.0 m', ...  
    'Mass', '27 000 tonnes'};
displayVehicleData('Semisubmersible: 2 Pontoons and 4 Legs', ...
    vesselData, 'semisub.jpg', 3);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Semisubmersible');
    disp('Two pontoons and four legs rig')
    disp('DP system: MIMO nonlinear PID controller');
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end


