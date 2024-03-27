echo on
% SIMcontainer - User editable script for simulating the dynamics of a 
% container ship under feedback control. The script concurrently simulates 
% the ship using both a linear model, defined in 'Lcontainer.m', and a 
% nonlinear model, defined in 'container.m'. The outcomes of both 
% simulations are then plotted side by side for comparative analysis.
%
% Calls:       container.m
%              Lcontainer.m
%              euler2.m
%
% Author:      Thor I. Fossen
% Date:        2018-07-21
% Revisions: 

echo off
disp('Simulating container.m and Lcontainer.m under PD control with psi_ref = 5 deg ...')

t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)

Kp = 1;      % controller P gain
Td = 10;     % controller derivative time

% initial states:
x1 = [7 0 0 0 0 0 0 0 0 70]';   % x = [ u v r x y psi p phi delta n ]'
x2 = [7 0 0 0 0 0 0 0 0]';     

%% MAIN LOOP
N = round(t_f/h);                       % number of samples
simdata1 = zeros(N+1,length(x1)+2);     % memory allocation
simdata2 = zeros(N+1,length(x2)+2);     % memory allocation
                
for i=1:N+1

    time = (i-1) * h;                   % simulation time in seconds

    r   = x1(3);
    psi = x1(6);
    
    % control system
    psi_ref = deg2rad(5);                           % desired heading
    delta_c = -Kp * ( ssa(psi-psi_ref) + Td * r);   % PD controller
    n_c = 70;
    
    % ship model
    [xdot1,U1] = container(x1,[delta_c n_c]);       % ship models
    [xdot2,U2] = Lcontainer(x2,delta_c);
   
    % store data for presentation
    simdata1(i,:) = [time,x1',U1]; 
    simdata2(i,:) = [time,x2',U2]; 
    
    % numerical integration
    x1 = euler2(xdot1,x1,h);                         % Euler's method
    x2 = euler2(xdot2,x2,h);   

end

%% PLOTS
t1     = simdata1(:,1);
u1     = simdata1(:,2); 
v1     = simdata(:,3);          
r1     = rad2deg(simdata1(:,4));   
x1     = simdata1(:,5);
y1     = simdata1(:,6);
psi1   = rad2deg(simdata1(:,7));
p1     = rad2deg(simdata1(:,8));
phi1   = rad2deg(simdata1(:,9));
delta1 = rad2deg(simdata1(:,10));
n1     = simdata1(:,11);
U1     = simdata1(:,12);

t2     = simdata2(:,1);
u2     = simdata2(:,2); 
v2     = simdata2(:,3);          
r2     = rad2deg(simdata2(:,4));   
x2     = simdata2(:,5);
y2     = simdata2(:,6);
psi2   = rad2deg(simdata2(:,7));
p2     = rad2deg(simdata2(:,8));
phi2   = rad2deg(simdata2(:,9));
delta2 = rad2deg(simdata2(:,10));
U2     = simdata2(:,11);

figure(1)
plot(y1,x1,'r',y2,x2,'b')
grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position')
legend('nonlinear model','linear model','Location','best')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2)
subplot(221),plot(t1,r1,'r',t,r2,'b'),xlabel('time (s)'),title('yaw rate r (deg/s)'),grid
legend('nonlinear model','linear model','Location','best')
subplot(222),plot(t1,phi1,'r',t,phi2,'b'),xlabel('time (s)'),title('roll angle \phi (deg/s)'),grid
legend('nonlinear model','linear model','Location','best')
subplot(223),plot(t1,psi1,'r',t,psi2,'b'),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
legend('nonlinear model','linear model','Location','best')
subplot(224),plot(t1,delta1,'r',t,delta2,'b'),xlabel('time (s)'),title('rudder angle \delta (deg)'),grid
legend('nonlinear model','linear model','Location','best')
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)