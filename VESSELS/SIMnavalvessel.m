function SIMnavalvessel()
% SIMnavalvessel is compatibel with MATLAB and GNU Octave (www.octave.org). 
% Thisscript simulates a naval vessel under PD heading control.
%
% Dependencies:  
%   navalvessel.m   - Naval vessel dynamics (Blanke and Christensen, 1993).
%   euler2.m        - Euler's integration method.
%
% Reference: 
%   M. Blanke and A. Christensen (1993). Rudder-roll damping autopilot 
%   robustness to sway-yaw-roll couplings. 10th Ship Control 
%   Systems Symposium, Ottawa, Canada.
%
% Author:      Thor I. Fossen
% Date:        2019-05-27
% Revisions:
%   2024-03-27 : Added animation of the ship North-East positions.
%   2024-04-22 : Enhanced compatibility with GNU Octave.

clear animateShip   % clear the persistent animation variables

t_f = 600;          % final simulation time (sec)
h   = 0.05;         % sample time (sec)
N = round(t_f/h);   % number of samples

Kp = 1e6;           % controller proportional gain
Td = 20;            % controller derivative time

% Initial states
x  = [6 0 0 0 0 0 ]';   % x = [u v p r phi psi] ]'   
pos = [0 0 0]';         % initial position expressed in NED

%% MAIN LOOP
simdata = zeros(N+1,11);                % memory allocation

for i=1:N+1
    time = (i-1)*h;                     % simulation time in seconds

    % Measurements
    r   = x(4) + 0.001 * randn;
    psi = x(6) + 0.01 * randn;
    
    % Control system
    psi_ref = deg2rad(20);                          % desired heading
    tauX = 1e5;                                    % surge command 
    tauN = -Kp * ( ssa(psi - psi_ref) + Td * r );  % PD heading controller
    
    % Ship dynamics
    tau = [tauX 0 0 tauN]';
    [xdot,U] = navalvessel(x,tau);       
   
    % Store data for presentation
    simdata(i,:) = [ time, x', tauN, U, pos(1:2)' ]; 
    
    % Numerical integration
    x = euler2(xdot,x,h);                 
    pos = pos + h * Rzyx(0,0,psi) * [x(1) x(2) x(4)]'; 

end

%% PLOTS
scSz = get(0, 'ScreenSize'); % Returns [left bottom width height]

t     = simdata(:,1);
u     = simdata(:,2); 
v     = simdata(:,3);          
p     = rad2deg(simdata(:,4));   
r     = rad2deg(simdata(:,5));
phi   = rad2deg(simdata(:,6));
psi   = rad2deg(simdata(:,7));
tauN  = simdata(:,8);
U     = simdata(:,9);
x     = simdata(:,10);
y     = simdata(:,11);

% Plot and animation of the North-East positions
figure(1)
if isoctave() % Octave NE-plot
    plot(y,x,'b')
    xlabel('East'); ylabel('North');title('North-East plot (m)')
    grid,axis('equal')
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
else % Matlab animation
    shipSize = 0.2;
    set(gcf, 'Position', [1, 1, 0.5*scSz(3), scSz(4)]);
    animateShip(x,y,shipSize,'b-',1);
end

figure(2)
plot(t,U,'b')
xlabel('time (s)'),title('Speed U (m/s)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(3)
subplot(221)
plot(t,r,'b'),xlabel('time (s)'),title('Yaw rate r (deg/s)'),grid
subplot(222)
plot(t,phi,'b'),xlabel('time (s)'),title('Roll angle \phi (deg)'),grid
subplot(223)
plot(t,psi,'b',[0,t(end)],rad2deg([psi_ref psi_ref])),xlabel('time (s)'),
title('Yaw angle \psi (deg)'),grid
subplot(224)
plot(t,tauN,'b'),xlabel('time (s)'),title('Yaw control moment (Nm)'),grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

end
