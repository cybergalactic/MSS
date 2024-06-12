% ExHybridPath is compatibel with MATLAB and GNU Octave (www.octave.org).
% The function calculates a hybrid continuous path based through specified
% waypoints. Use the example file 'ExHybridPath.m' to generate and test
% the hybrid path.
%
% Dependencies:
%   getPathSignal.m - Generates the coefficients for subpaths between 
%                     given waypoints
%
% Author:        Roger Skjetne
% Date created:  2014-02-10  Roger Skjetne.
% Revised:       
%   2024-06-12 : Added compability to GNU Octave. 

clearvars;

%% Initialization
X_wp = [0 5 10 15 20 25 30 35 40];      % Waypoints
Y_wp = [0 9 10 11 20 10 20 10 20];
% X_wp = [10  5  0  5 10  5  0  5 10  5  0  5 10];
% Y_wp = [10 15 10  5 10 15 10  5 10 15 10  5 10];

s           = 1.0;  % Path parameter in [0, NumSubPaths]
lambda      = 0.25; % Curvature constant.
r           = 2;    % Differentiability order.
PlotHandle  = 1;    % Empty if no plotting

%% Run scripts
WP = [X_wp' Y_wp'];
Path = hybridPath(WP,r,lambda,PlotHandle);
pathSignals = getPathSignals(Path,s);

%% Plotting
if PlotHandle
    figure(PlotHandle); hold on;
    plot(pathSignals.pd(1),pathSignals.pd(2),'sm','LineWidth',1.5);
    plot([pathSignals.pd(1) pathSignals.pd(1)+pathSignals.pd_der{1}(1)],...
        [pathSignals.pd(2) pathSignals.pd(2)+pathSignals.pd_der{1}(2)],'m','LineWidth',1.25);
    s = 0:0.01:Path.NumSubpaths;

    for k=1:Path.Order
        pd_der = zeros(length(s),2);
        for j=1:length(s)
            PS = getPathSignals(Path,s(j));
            pd_der(j,:) = PS.pd_der{k};
        end

        % Plotting derivatives
        figure(PlotHandle+k); clf;
        plot(s,pd_der(:,1),s,pd_der(:,2),'LineWidth',1.25);
        text_x = ['x_d^{s^',num2str(k),'}'];
        text_y = ['y_d^{s^',num2str(k),'}'];
        legend(text_x,text_y); grid on;
    end

end
