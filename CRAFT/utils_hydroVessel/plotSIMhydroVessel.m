function plotSIMhydroVessel(t, simdata, environment)

legendLocation = 'best';
if isoctave; legendLocation = 'northeast'; end

xn    = simdata(:,1);
yn    = simdata(:,2);
zn    = simdata(:,3);
phi   = ssa(simdata(:,4));
theta = ssa(simdata(:,5));
psi   = ssa(simdata(:,6));

u     = simdata(:,7);
v     = simdata(:,8);
w     = simdata(:,9);
p     = simdata(:,10);
q     = simdata(:,11);
r     = simdata(:,12);

tau           = simdata(:,13:18);
tau_wave1     = simdata(:,19:24);
waveElevation = simdata(:,25);

spreadingFlag = environment.spreadingFlag;
spectrumType = environment.spectrumType;
Hs = environment.Hs;
w0 = environment.w0;
beta_wave = environment.beta_wave;
mu = environment.mu;
Omega = environment.Omega;
S_M = environment.S_M;

U = sqrt( u.^2 + v.^2 );        % Vessel speed (m/s)

%% Position and Euler angle plots
figure(2); clf;
figure(gcf)
subplot(321),plot(yn,xn)
xlabel('East (m)')
ylabel('North (m)')
title('North-East positions (m)'),grid
subplot(322),plot(t,zn)
xlabel('time (s)'),title('Down position (m)'),grid
subplot(312),plot(t,rad2deg(phi),t,rad2deg(theta))
xlabel('time (s)'),title('Roll and pitch angles (deg)'),grid
legend('Roll angle (deg)','Pitch angle (deg)')
subplot(313),plot(t,rad2deg(psi))
xlabel('time (s)'),title('Heading angle (deg)'),grid
legend('Yaw angle (deg)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

%% Velocity plots
figure(3); clf;
figure(gcf)
subplot(311),plot(t,U)
xlabel('time (s)'),title('Speed (m/s)'),grid
subplot(312),plot(t,u,t,v,t,w)
xlabel('time (s)'),title('Linear velocities (m/s)'),grid
legend('u (m/s)','v (m/s)','w (m/s)')
subplot(313),plot(t,rad2deg(p),t,rad2deg(q),t,rad2deg(r))
xlabel('time (s)'),title('Angular velocities (deg/s)'),grid
legend('p (deg/s)','q (deg/s)','r (deg/s)','Location',legendLocation)
set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

%% Plot the 6-DOF control forces
figure(4); clf;
figure(gcf)
DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)',...
    'Roll (Nm)', 'Pitch (Nm)', 'Yaw (Nm)'};
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, tau(:, DOF), 'LineWidth', 1.5);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF});
end

if ~isoctave
    sgtitle(['Generalized control forces'],'FontSize', 12);
end

%% Plot the 6-DOF 1st-order wave forces
figure(5); clf;
figure(gcf)
DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)',...
    'Roll (Nm)', 'Pitch (Nm)', 'Yaw (Nm)'};
for DOF = 1:6
    subplot(6, 1, DOF);
    plot(t, tau_wave1(:, DOF), 'LineWidth', 1.5);
    xlabel('Time (s)');
    grid on;
    legend(DOF_txt{DOF});
end

if ~isoctave
    sgtitle(['Generalized 1st-order wave-induced forces for \beta_{wave} = ' ...
        num2str(rad2deg(beta_wave)), '°, H_s = ', num2str(Hs), ' m and ' ...
        '\omega_0 = ', num2str(w0), ' rad/s'],'FontSize', 12);
end

%% Plot the wave spectrum and wave elevation
figure(6); clf;
figure(gcf)
subplot(211);
hold on;

if spreadingFlag
    % Plot the wave spectrum for the specific directions
    hold on;
    plot(Omega, S_M(:, floor(length(mu)/2)),'b->','Markersize',5,'LineWidth',1.5);
    plot(Omega, S_M(:, floor(length(mu)/4)),'k-o','Markersize',5, 'LineWidth',1.5);
    plot(Omega, S_M(:, length(mu)),'g-', 'LineWidth',2.0);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))],'r-.', 'LineWidth', 1.5)
    legend('\mu = 0 deg', '\mu = 45 deg', '\mu = 90 deg',...
        ['\omega_0 = ', num2str(w0), ' rad/s']);
    hold off;
else
    hold on
    plot(Omega, S_M(:, 1), 'b-', 'LineWidth', 1.5);
    plot([w0, w0], [min(min(S_M)), max(max(S_M))], 'r-.','LineWidth', 1.5)
    legend('S(\omega)', ['\omega_0 = ', num2str(w0), ' rad/s']);
    hold off
end

xlabel('Omega (rad/s)');
ylabel('m^2 s');
title([spectrumType, ' spectrum']);
grid on;

% Plot the wave elevation
subplot(212);
plot(t, waveElevation, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('m');
title(['Wave elevation for wave direction \beta_{wave} = ', ...
    num2str(rad2deg(beta_wave)), '°, H_s = ', num2str(Hs), ' m and ' ...
    '\omega_0 = ', num2str(w0), ' rad/s']);
grid on;

set(findall(gcf,'type','line'),'linewidth',1.5)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',12)

end