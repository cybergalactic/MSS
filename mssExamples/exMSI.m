% exMSI is compatible with MATLAB and GNU Octave (www.octave.org).
% The script plots the ISO 2631-1 (1997) and O'Hanlon and McCauley (1974)
% Motion Sickness Incidence curves
%
% Author:    Thor I. Fossen
% Date:      2001-11-05
% Revisions: 

%% Plot the O'Hanlon and McCauley MSI
w_0 = 0.01:0.1:2;       % Wave frequency (rad/s)
U = 15;                 % Vessel speed (m/s)
beta = deg2rad(180);    % Wave encounter angle (rad)

a_z = [0.5 1 2 3 4 5];  % Vertical accelerations (m/s^2)
idx = ['x','o','d','h','s','v'];

figure(1); clf;
hold on

for i = 1:6
    w_e = encounter(beta, U, w_0); % Frequency of encounter (rad/s)
    msi = HMmsi(a_z(i), w_e);      % O'Hanlon and McCauley MSI (%)
    plot(w_e, msi, ['-', idx(i)]); % Plot curve
end

hold off
title('O''Hanlon and McCauley (1974) Motion Sickness Incidence (MSI)')
xlabel('Frequency of encounter \omega_e (rad/s)')
ylabel('MSI (%)')
legend('a_z = 0.5 (m/s^2)', 'a_z = 1 (m/s^2)', 'a_z = 2 (m/s^2)', ...
    'a_z = 3 (m/s^2)', 'a_z = 4 (m/s^2)', 'a_z = 5 (m/s^2)');
title('Percentage of Persons that become Seasick during a 2 Hours Voyage', ...
    'O''Hanlon and McCauley (1974) Motion Sickness Incidence (MSI)');
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

%% Plot the ISO MSI
clearvars

t = [0.5 1 2 4 8];
idx = ['^','*','d','s','o'];
figure(2); clf;
hold on

for i = 1:5
    [a_z, w_e] = ISOmsi(t(i));     % ISO 2631 MSI
    plot(w_e, a_z,['-', idx(i)]);  % Plot curve
end

hold off
title('ISO 2631 Motion Sickness Incidence (MSI)')
xlabel('Frequency of encounter \omega_e (rad/s)')
ylabel('a_z (m/s^2)')
legend('ISO 30 MIN', 'ISO 1 HR', 'ISO 2 HRS', 'ISO 4 HRS', 'ISO 8 HRS')
grid
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)