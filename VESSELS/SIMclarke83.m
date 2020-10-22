% SIMclarke83 User editable script for simulation of a ship with main
% dimensions L, B and T. The hydrodynamic data is computed using clarke83.m.
%
% Author:      Thor I. Fossen
% Date:        2020-10-22
% Revisions: 

clear all

h = 0.05;               % sampling time
N = 10000;              % number of samples

psi_ref = 10 * pi/180;  % Heading angle setpoint
w_n = 0.1;              % Bandwidth
Kp = w_n^2;             % Proportional gain 
Kd = 2 * w_n;           % Derivative  gain 

% initial values
eta = zeros(3,1);       % x, y, psi
nu  = [0 0 0 ]';        % u, v, r

% ship model
L = 100;    % length (m)
B = 20;     % beam (m)
T = 10;     % draft (m)
Cb = 0.8;   % block coefficient, Cb = V / (L*B*T) where V is the displaced volume
R66 = 0.27*L; % radius of gyration (smaller vessels R66 ≈ 0.25L, tankers R66 ≈ 0.27L)
xg = -3;    % x-coordinate of the CG   
     
% *************** MAIN SIMULATION LOOP ************************
for i=1:N+1
    t = (i-1)*h;                   % simulation time in seconds
  
    % maneuvering model
    [M,C,D] = clarke83(nu,L,B,T,Cb,R66,xg);
    
    % Control system (constant thrust + PD heading controller)
    tau = [1000000
           0
           M(3,3) * ( Kp * ssa(psi_ref-eta(3)) - Kd * nu(3) ) ]; 
     
    % Differential equations
    etadot = Rzyx(0,0,eta(3)) * nu;
    nudot = inv(M) * (tau - C * nu - D * nu);
    
    % store data for presentation
    xout(i,:) = [t, eta', nu'];
    
    % numerical integration
    eta = euler2(etadot, eta, h);             % Euler integration
    nu  = euler2(nudot, nu, h);     

end

% *************** END SIMULATION LOOP ************************

% Time-series
t     = xout(:,1);
x     = xout(:,2); 
y     = xout(:,3);            
psi   = xout(:,4) * 180/pi;
u     = xout(:,5);
v     = xout(:,6);
r     = xout(:,7) * 180/pi; 
U     = sqrt(u.^2 + v.^2);

% plots
figure(1)
subplot(311)
plot(y,x); grid,axis('equal')
xlabel('East'),ylabel('North'),title('Ship position (m)')

subplot(312)
plot(t,U)
xlabel('time (s)'),title('Ship speed (m/s)'),grid

subplot(313)
plot(t,psi)
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid



     
