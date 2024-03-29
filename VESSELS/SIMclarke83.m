% SIMclarke83 User editable script for simulation of a ship with main
% dimensions L, B, and T. The hydrodynamic data is based on:
%
% Reference: CLARKE, D., GEDLING, P. and HINE. G. (1983). The application of 
% manoeuvring criteria in hull design using linear thory. Trans.  
% R. lnsm nav. Archit.  125, 45-68. 
%
% Calls:     clarke83.m
%
% Author:    Thor I. Fossen
% Date:      2020-10-22
%            2024-03-27 Using forward and backward Euler to integrate xdot
%                       Added animation of the ship North-East positions

clear animateShip  % clear the persistent animation variables
clearvars;

h = 0.05;               % sampling time
N = 10000;              % number of samples

psi_ref = deg2rad(10);  % heading angle setpoint
w_n = 0.1;              % closed-loop natural frequency
Kp = w_n^2;             % proportional gain 
Kd = 2 * w_n;           % derivative  gain 

% Initial values
eta = zeros(3,1);       % x, y, psi
nu  = [0 0 0 ]';        % u, v, r

% Ship model
L = 100;      % length (m)
B = 20;       % beam (m)
T = 10;       % draft (m)
Cb = 0.8;     % block coefficient, Cb = V / (L*B*T) where V is the displaced volume
R66 = 0.27*L; % radius of gyration (smaller vessels R66 ≈ 0.25L, tankers R66 ≈ 0.27L)
xg = -3;      % x-coordinate of the CG   
     
% MAIN LOOP
simdata = zeros(N+1,7);                   % table for simulation data

for i=1:N+1

    t = (i-1)*h;                   % simulation time in seconds

    % Linear maneuvering model
    U = sqrt(nu(1)^2 + nu(2)^2);
    [M,N] = clarke83(U,L,B,T,Cb,R66,xg);

    % Control system (constant thrust + PD heading controller)
    tau = [1000000
        0
        M(3,3) * ( Kp * ssa(psi_ref-eta(3)) - Kd * nu(3) ) ];

    % Differential equations
    etadot = Rzyx(0,0,eta(3)) * nu;
    nudot = M \ (tau - N * nu);

    % Store data for presentation
    simdata(i,:) = [t, eta', nu'];

    % Euler's integration methods (k+1), (Fossen 2021, Eq. B27-B28)
    % x = x + h * xdot is replaced by forward and backward Euler integration
    nu  = nu + h * nudot;                       % Forward Euler
    eta = eta + h * Rzyx(0,0,eta(3)) * nu;      % Backward Euler

end

%% PLOTS
t     = simdata(:,1);
x     = simdata(:,2); 
y     = simdata(:,3);            
psi   = rad2deg(simdata(:,4));
u     = simdata(:,5);
v     = simdata(:,6);
r     = rad2deg(simdata(:,7)); 
U     = sqrt(u.^2 + v.^2);

figNo = 1;
shipSize = 0.6;
animateShip(x,y,shipSize,'b-',figNo); % animation of the North-East positions

figure(2)
subplot(211)
plot(t,U)
xlabel('time (s)'),title('Ship speed (m/s)'),grid
subplot(212)
plot(t,psi)
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)


     
