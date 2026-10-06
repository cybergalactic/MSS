function plotBv(vessel)
% plotBv plots the zero-speed potential damping and the constant
%    power-based matrices B_eq, Bv and D = B_eq + Bv.
%
%    plotBv(vessel)
%
% Input: 
%    vessel:  MSS vessel structure 
%
% Author:    Thor I. Fossen
% Date:      2020-03-08 
% Revisions: 2026-09-26 Added constant power-based damping overlays.
%            2026-09-27 Removed the separate seakeeping vessel.Bv model.
%                       Plot B(inf) as a separate red asterisk.

w       = vessel.freqs;
Nfreq   = length(w);
velno   = 1;
B       = vessel.B;

% VERES uses omega = 10 rad/s to represent infinite frequency. Do not
% connect this artificial point to the physical frequency-domain curves.
infIdx  = abs(w - 10) < 10*eps(10);
lineIdx = ~infIdx;

% Constant power-based damping matrices are plotted from 0 to 10 rad/s
wConst = [0; 10];

% Check that the power-based damping matrices are available
hasPowerBeq = isfield(vessel, 'powerBased') && ...
    isstruct(vessel.powerBased) && ...
    isfield(vessel.powerBased, 'B_eq') && ...
    ~isempty(vessel.powerBased.B_eq);

hasPowerBv = isfield(vessel, 'powerBased') && ...
    isstruct(vessel.powerBased) && ...
    isfield(vessel.powerBased, 'Bv') && ...
    ~isempty(vessel.powerBased.Bv);

if ~hasPowerBeq || ~hasPowerBv
    error(['Run computeManeuveringModel before plotBv so that ', ...
        'vessel.powerBased.B_eq and vessel.powerBased.Bv are available.']);
end

B_eq = vessel.powerBased.B_eq(:,:,1,1);
if size(B_eq,1) ~= 6 || size(B_eq,2) ~= 6
    error('The power-based B_eq matrix must have size 6-by-6.');
end

Bv_eq = vessel.powerBased.Bv(:,:,1,1);
if size(Bv_eq,1) ~= 6 || size(Bv_eq,2) ~= 6
    error('The power-based Bv matrix must have size 6-by-6.');
end

figno = 50;

% Longitudinal plots
k = 1;
figure(figno)

for i = 1:2:5
    for j = 1:2:5

        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        splot = 330+k;
        subplot(splot)

        % Potential damping at physical frequencies
        hB = plot(w(lineIdx),Bplot(lineIdx),'b-o');
        hold on
        plotHandles = hB;
        plotNames = {'B (potential)'};

        % Infinite-frequency potential damping stored at omega = 10 rad/s
        if any(infIdx)
            hInf = plot(w(infIdx),Bplot(infIdx),'ro', ...
                'MarkerFaceColor','b','MarkerSize',6);
            plotHandles(end+1) = hInf;
            plotNames{end+1} = 'B(\infty)';
        end

        % Power-based equivalent potential damping
        hBeq = plot(wConst,B_eq(i,j)*ones(2,1), 'c-.','linewidth',2);
        plotHandles(end+1) = hBeq;
        plotNames{end+1} = 'B_{eq} (power-based)';

        % Power-based viscous damping
        hBv = plot(wConst,Bv_eq(i,j)*ones(2,1), 'r-','linewidth',2);
        plotHandles(end+1) = hBv;
        plotNames{end+1} = 'B_v (power-based)';

        % Total power-based damping
        D_eq = B_eq(i,j) + Bv_eq(i,j);
        hDeq = plot(wConst,D_eq*ones(2,1), 'm-.','linewidth',2);
        plotHandles(end+1) = hDeq;
        plotNames{end+1} = 'B_{eq}+B_v (power-based)';

        hold off
        grid
        legend(plotHandles,plotNames,'Location','best','FontSize',8)
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');
        xlabel('frequency (rad/s)')
        title(Hw);

        k = k + 1;

    end
end

% Lateral plots
figno = figno + 1;
k = 1;
figure(figno)

for i = 2:2:6
    for j = 2:2:6

        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        splot = 330+k;
        subplot(splot)

        % Potential damping at physical frequencies
        hB = plot(w(lineIdx),Bplot(lineIdx),'b-o');
        hold on
        plotHandles = hB;
        plotNames = {'B (potential)'};

        % Infinite-frequency potential damping stored at omega = 10 rad/s
        if any(infIdx)
            hInf = plot(w(infIdx),Bplot(infIdx),'ro', ...
                'MarkerFaceColor','b','MarkerSize',6);
            plotHandles(end+1) = hInf;
            plotNames{end+1} = 'B(\infty)';
        end

        % Power-based equivalent potential damping
        hBeq = plot(wConst,B_eq(i,j)*ones(2,1), 'c-.','linewidth',2);
        plotHandles(end+1) = hBeq;
        plotNames{end+1} = 'B_{eq} (power-based)';

        % Power-based viscous damping
        hBv = plot(wConst,Bv_eq(i,j)*ones(2,1), 'r-','linewidth',2);
        plotHandles(end+1) = hBv;
        plotNames{end+1} = 'B_v (power-based)';

        % Total power-based damping
        D_eq = B_eq(i,j) + Bv_eq(i,j);
        hDeq = plot(wConst,D_eq*ones(2,1), 'm-.','linewidth',2);
        plotHandles(end+1) = hDeq;
        plotNames{end+1} = 'B_{eq}+B_v (power-based)';

        hold off
        grid
        legend(plotHandles,plotNames,'Location','best','FontSize',8)
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');
        xlabel('frequency (rad/s)')
        title(Hw);

        k = k + 1;

    end
end

end