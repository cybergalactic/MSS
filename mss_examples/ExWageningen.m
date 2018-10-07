% ExWageningen computes thrust and torque curves for a propeller using 
%              the Wageningen B-series data. 
%
% The B-series propellers were designed and tested at the Netherlands Ship
% Model Basin in Wageningen. The open_water characteristics of 120
% propeller models of the B_series were tested at N.S.M.B. and analyzed
% with multiple polynomial regression analysis. The derived polynomials
% express the thrust and torque coefficients in terms of the number of
% blades, the blade area ratio, the pitch_diameter ratio and the advance
% coefficient.
%
% Reference: Barnitsas, M.M., Ray, D. and Kinley, P.
% Kt, Kq and Efficiency Curves for the Wageningen B-Series Propellers
% https://deepblue.lib.umich.edu/handle/2027.42/3557
%
% Author:    Thor I. Fossen
% Date:      7 October 2018
% Revisions:
clear all

rho = 1025;   % Density of water (kg/m^3)
D = 5;        % Propeller diameter (m)
PD = 1.4;     % pitch/diameter ratio (typically 0.5-2.5)
AEAO = 0.65;  % blade area ratio (ratio of expanded blade area to propeller disk area)
z = 4;        % number of propeller blades

% Comput K_T and K_Q for advance velocites J_a
J_a = -1:0.01:1;
for i = 1:length(J_a)
    [K_T(i), K_Q(i)] = wageningen(J_a(i),PD,AEAO,z);
end

% Compute K_T and K_Q for J_a = 0 (Bollard pull)
[K_T0, K_Q0] = wageningen(0,PD,AEAO,z);

% Compute thrust as a funtion of propeller 
n = 0:1:10;                           % propeller [RPS]
T = rho * D^4 * K_T0 * n .* abs(n);   % thrust [N]
    
%% PLots 
figure(gcf)
subplot(211)
plot(J_a,K_T,J_a,10*K_Q,'LineWidth',2)
grid on
xlabel('J_a - Advance Velocity');
legend('K_T','10 K_Q');
title(['Wageningen B-series propeller with ',num2str(z),' blades and P/D = ',num2str(PD), ...
    ', AE/AO = ', num2str(AEAO) ]);

subplot(212)
plot(n * 60 ,T/1000,'LineWidth',2)
grid on
xlabel('Propeller [RPM]');
ylabel('Thrust T [kN]');
title('Thrust in kN vs. proepeller RPM');

