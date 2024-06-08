function SIMclarke83()
% SIMclarke83 is compatible with MATLAB and GNU Octave (www.octave.org). 
% This script simulates a ship, characterizing its dynamics based on 
% specified main dimensions: length (L), breadth (B), and draft (T). It 
% uses hydrodynamic data based on:
%
% Reference:
%   D. Clarke, P. Gedling, and G. Hine (1983). The application of manoeuvring
%   criteria in hull design using linear theory. Transactions of the Royal
%   Institution of Naval Architects, Vol. 125, pp. 45-68.
%
% Dependencies:
%   clarke83.m - Function implementing Clarke's linear maneuvering model.
%
% Author: Thor I. Fossen
% Date:   2020-10-22
% Revisions: 
%   2024-03-27 : Using forward and backward Euler to integrate xdot.
%                Added animation of the ship North-East positions.
%   2024-04-19 : Enhanced compatibility with GNU Octave.

clear animateShip       % Clear the persistent animation variables
clearvars;
close all;

%% USER INPUTS
h = 0.05;               % Sampling time
N = 10000;              % Number of samples

psi_ref = deg2rad(10);  % Heading angle setpoint
w_n = 0.1;              % Closed-loop natural frequency
Kp = w_n^2;             % Proportional gain
Kd = 2 * w_n;           % Derivative  gain

% Initial values
eta = zeros(3,1);       % x, y, psi
nu  = [0 0 0 ]';        % u, v, r

% Ship model
L = 100;      % Length (m)
B = 20;       % Beam (m)
T = 10;       % Draft (m)
Cb = 0.8;     % Block coefficient, Cb = V / (L*B*T) where V is the displaced volume
R66 = 0.27*L; % Radius of gyration (smaller vessels R66 ≈ 0.25L, tankers R66 ≈ 0.27L)
xg = -3;      % x-coordinate of the CG

% Display simulation options
displayControlMethod();

%% MAIN LOOP
simdata = zeros(N+1,7);                   % table for simulation data

for i=1:N+1

    t = (i-1) * h;                        % simulation time in seconds

    % Linear maneuvering model
    U = sqrt(nu(1)^2 + nu(2)^2);
    [M,N] = clarke83(U,L,B,T,Cb,R66,xg);

    % Control system (constant thrust + PD heading controller)
    tau = [1000000
        0
        M(3,3) * ( Kp * ssa(psi_ref-eta(3)) - Kd * nu(3) ) ];

    % Differential equations
    nudot = M \ (tau - N * nu);

    % Store data for presentation
    simdata(i,:) = [t, eta', nu'];

    % Euler's integration methods (k+1), (Fossen 2021, Eq. B27-B28)
    % x = x + h * xdot is replaced by forward and backward Euler integration
    nu  = nu + h * nudot;                       % Forward Euler
    eta = eta + h * Rzyx(0,0,eta(3)) * nu;      % Backward Euler

end

%% PLOTS
scrSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

% Simdata(i,:) = [t, eta', nu']
t     = simdata(:,1);
x     = simdata(:,2);
y     = simdata(:,3);
psi   = rad2deg(simdata(:,4));
u     = simdata(:,5);
v     = simdata(:,6);
r     = rad2deg(simdata(:,7));
U     = sqrt(u.^2 + v.^2);

% Plot and animation of the North-East positions
figure(1)
shipSize = 0.5;
set(gcf, 'Position', [1, 1, 0.4*scrSz(3), scrSz(4)]);
animateShip(x,y,shipSize,'b-',1);

% Ship speed and yaw angle
figure(2);
subplot(211)
plot(t,U)
xlabel('Time (s)'),title('Ship speed (m/s)'),grid
subplot(212)
plot(t,psi)
xlabel('Time (s)'),title('Yaw angle \psi (deg)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)

% Display the vessel data and an image of the vessel
vesselData = {...
    'Length', [num2str(L) ' m'], ...
    'Beam', [num2str(B) ' m'], ...
    'Draft', [num2str(T) ' m'],... 
    'Block coefficient', num2str(Cb),...    
    'Mass', [num2str(M(1,1)/1000) ' tonnes']};
displayVehicleData('Ship characterized by L, B and T',...
    vesselData, 'clarke.jpg', 3);

end

%% DISPLAY CONTROL METHOD
function displayControlMethod()
    disp('--------------------------------------------------------------------');
    disp('MSS toolbox: Ship characterized by length (L), beam (B) and draft (T)');
    disp('Heading autopilot: PD control law')
    disp('Speed: Constant thrust')
    disp('--------------------------------------------------------------------');
    disp('Simulating...');
end

