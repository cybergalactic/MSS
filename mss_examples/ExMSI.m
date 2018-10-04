% ExMSI Plots the ISO 2631-1 (1997) and O'Hanlon and McCauley (1974)
%       Motion Sickness Incidence curves
%
% Author:    Thor I. Fossen
% Date:      5th November 2001
% Revisions: 

clf
conversion       % load conversion factors

%---------------------------------------------------------------------------
% Plot O'Hanlon and McCauley MSI
%---------------------------------------------------------------------------

w_0 = 0.01:0.1:2;  % wave frequency (rad/s)
U = 15;            % vessel speed (m/s)
beta = 180*D2R;    % wave encounter angle (rad)

a_z = [0.5 1 2 3 4 5]; % vertical accelerations (m/s^2)
idx = ['x','o','d','h','s','v'];

figure(1); hold on
for i = 1:6,
    w_e = encounter(w_0,U,beta);           % frequency of encounter (rad/s)
    msi = HMmsi(a_z(i),w_e);               % O'Hanlon and McCauley MSI (%)
    h=plot(w_e,msi,'b',w_e,msi,idx(i));    % plot curve
    hh(i)=h(2);                            % figure handle
end

hold off
title('O''Hanlon and McCauley (1974) Motion Sickness Incidence (MSI)')
xlabel('Frequency of encounter \omega_e (rad/s)')
ylabel('MSI (%)')
legend(hh,'a_z = 0.5 (m/s^2)','a_z  = 1 (m/s^2)','a_z  = 2 (m/s^2)','a_z  = 3 (m/s^2)',...
    'a_z  = 4 (m/s^2)','a_z  = 5 (m/s^2)');
grid

%---------------------------------------------------------------------------
% Plot ISO MSI
%---------------------------------------------------------------------------
clear all


t = [0.5 1 2 4 8];
idx = ['^','*','d','s','o'];
figure(2); hold on

for i = 1:5,
    [a_z,w_e] = ISOmsi(t(i));              % ISO 2631 MSI
    h = plot(w_e,a_z,'b',w_e,a_z,idx(i));  % plot curve
    hh(i)=h(2);                            % figure handle
end

hold off
title('ISO 2631 Motion Sickness Incidence (MSI)')
xlabel('Frequency of encounter \omega_e (rad/s)')
ylabel('a_z (m/s^2)')
legend(hh,'ISO 30 MIN','ISO 1 HR','ISO 2 HRS','ISO 4 HRS','ISO 8 HRS')
grid

