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
% Reference: Barnitsas, M.M., Ray, D. and Kinley, P. (1981).
% KT, KQ and Efficiency Curves for the Wageningen B-Series Propellers
% http://deepblue.lib.umich.edu/handle/2027.42/3557
%
% Author:    Thor I. Fossen
% Date:      7 October 2018
% Revisions: 18 June 2019 - added curve fitting of KT and KQ data
clear all

rho = 1025;   % Density of water (kg/m^3)
D = 5;        % Propeller diameter (m)
PD = 1.4;     % pitch/diameter ratio (typically 0.5-2.5)
AEAO = 0.65;  % blade area ratio (ratio of expanded blade area to propeller disk area)
z = 4;        % number of propeller blades

% Comput KT and KQ for advance velocites Ja
Ja = -0.5:0.01:1.8;
for i = 1:length(Ja)
    [KT(i), KQ(i)] = wageningen(Ja(i),PD,AEAO,z);
    eta_0 = Ja(i)/(2*pi) * (KT(i)/KQ(i));
    if (KT(i) > 0 && KQ(i) >0 && Ja(i) >0)
        eta(i) = eta_0;
    else
        eta(i) = 0;
    end
end

    
% Compute KT and KQ for Ja = 0 (Bollard pull)
[KT_0, KQ_0] = wageningen(0,PD,AEAO,z);

% Compute thrust [N]
n = 0:0.1:10;                          % propeller [RPS]
T = rho * D^4 * KT_0 * n .* abs(n);    % thrust [N]

% Fit KT and KQ data to straight lines 
Jdata = 0:0.01:1.6;
for i = 1:length(Jdata)
    [KTdata(i), KQdata(i)] = wageningen(Jdata(i),PD,AEAO,z);
end

alpha = polyfit(Jdata,KTdata,1)    % KT = alpha(1)*Ja + alpha(2)
beta  = polyfit(Jdata,KQdata,1)    % KQ = beta(1) *Ja + beta(2)

%% PLots 
figure(gcf)
subplot(211)
plot(Ja,KT,'-k',Ja,10*KQ,'-.r',Ja(eta>0),eta(eta>0),':b','LineWidth',2)
hold on
plot(Jdata,alpha(1)*Jdata+alpha(2),'-k','LineWidth',2)
plot(Jdata,10*(beta(1)*Jdata+beta(2)),'-r','LineWidth',2)
hold off
xlabel('J_a - Advance velocity');
legend('K_T','10 K_Q','\eta_0'); grid on;
title(['Wageningen B-series propeller with ',num2str(z), ...
    ' blades and P/D = ',num2str(PD),', AE/AO = ', num2str(AEAO) ]);

subplot(212)
plot(n * 60 ,T/1000,'k','LineWidth',2); grid on;
xlabel('Propeller [RPM]');
ylabel('Thrust T [kN]');
title('Thrust in kN vs. propeller RPM');

