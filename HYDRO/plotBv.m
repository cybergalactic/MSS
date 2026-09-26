function plotBv(vessel)
% plotBv plots the zero-speed potential, viscous and total damping 
%    B_total(w) = B(w) + Bv(w) as a function of the frequency: 
%
%    plotBv(vessel)
%
% Input: 
%    vessel:  MSS vessel structure 
%
% Author:    Thor I. Fossen
% Date:      2020-03-08 First version
% Revisions: 2026-09-26 Added constant power-based damping overlays

w       = vessel.freqs;
Nfreq   = length(w);
velno   = 1;
B       = vessel.B;

% VERES uses omega = 10 rad/s to represent infinite frequency. Do not
% connect this artificial point to the physical frequency-domain curves.
infIdx = abs(w - 10) < 10*eps(10);
lineIdx = ~infIdx;

% Frequency-dependent viscous damping used by the Cummins model
Bv = vessel.Bv;

% Expand a constant Cummins Bv over all frequencies when necessary.
if size(Bv,1) ~= 6 || size(Bv,2) ~= 6
    error('The viscous damping matrix Bv must have size 6-by-6-by-Nfreq.');
elseif size(Bv,3) == 1
    Bv = repmat(Bv,1,1,Nfreq);
elseif size(Bv,3) ~= Nfreq
    error('The frequency-dependent Bv must contain one matrix per frequency.');
end

% Constant equivalent matrices used by the power-based model. These are
% plotted separately and are never combined with the Cummins matrices.
hasPowerBeq = isfield(vessel, 'powerBased') && ...
    isstruct(vessel.powerBased) && ...
    isfield(vessel.powerBased, 'B_eq') && ...
    ~isempty(vessel.powerBased.B_eq);

hasPowerBv = isfield(vessel, 'powerBased') && ...
    isstruct(vessel.powerBased) && ...
    isfield(vessel.powerBased, 'Bv') && ...
    ~isempty(vessel.powerBased.Bv);

if hasPowerBeq
    B_eq = vessel.powerBased.B_eq(:,:,1,1);
    if size(B_eq,1) ~= 6 || size(B_eq,2) ~= 6
        error('The power-based B_eq matrix must have size 6-by-6.');
    end
end

if hasPowerBv
    Bv_eq = vessel.powerBased.Bv(:,:,1,1);
    if size(Bv_eq,1) ~= 6 || size(Bv_eq,2) ~= 6
        error('The power-based Bv matrix must have size 6-by-6.');
    end
end

figno = 50;

% Longitudinal plots
k = 1;
figure(figno)
        
for i = 1:2:5
    for j = 1:2:5
        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        Bvplot = reshape(Bv(i,j,:),Nfreq,1);
        Btotalplot = Bplot + Bvplot;
        splot = 330+k;
        subplot(splot)
        hTotal = plot(w(lineIdx),Btotalplot(lineIdx),'k-','linewidth',2);
        hold on
        hB = plot(w(lineIdx),Bplot(lineIdx),'b-o');
        hBv = plot(w(lineIdx),Bvplot(lineIdx),'r-','linewidth',2);
        plotHandles = [hTotal hB hBv];
        plotNames = {'B+B_v (seakeeping)', 'B (potential)', ...
            'B_v(\omega) (viscous)'};

        if any(infIdx)
            hInf = plot(w(infIdx),Btotalplot(infIdx),'kx','linewidth',2);
            plot(w(infIdx),Bplot(infIdx),'bx','linewidth',2)
            plot(w(infIdx),Bvplot(infIdx),'rx','linewidth',2)
            plotHandles(end+1) = hInf;
            plotNames{end+1} = '\omega=\infty (stored at 10 rad/s)';
        end

        if hasPowerBeq
            Beqplot = B_eq(i,j)*ones(Nfreq,1);
            hBeq = plot(w(lineIdx),Beqplot(lineIdx),'c-.','linewidth',2);
            plotHandles(end+1) = hBeq;
            plotNames{end+1} = 'B_{eq} (power-based)';
            if any(infIdx)
                plot(w(infIdx),Beqplot(infIdx),'cx','linewidth',2)
            end
        end

        if hasPowerBeq && hasPowerBv
            D_eq = B_eq(i,j) + Bv_eq(i,j);
            Deqplot = D_eq*ones(Nfreq,1);
            hDeq = plot(w(lineIdx),Deqplot(lineIdx),'m-.','linewidth',2);
            plotHandles(end+1) = hDeq;
            plotNames{end+1} = 'B_{eq}+B_v (power-based)';
            if any(infIdx)
                plot(w(infIdx),Deqplot(infIdx),'mx','linewidth',2)
            end
        end

        hold off
        grid
        legend(plotHandles,plotNames,'Location','best','FontSize',8)
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');
        xlabel('frequency (rad/s)')
        title(Hw);
        k = k +1;
        
    end
    
end

% Lateral plots
figno = figno + 1;
k = 1;
figure(figno)

for i = 2:2:6  
    for j = 2:2:6       
        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        Bvplot = reshape(Bv(i,j,:),Nfreq,1);
        Btotalplot = Bplot + Bvplot;
        splot = 330+k;
        subplot(splot)
        hTotal = plot(w(lineIdx),Btotalplot(lineIdx),'k-','linewidth',2);
        hold on
        hB = plot(w(lineIdx),Bplot(lineIdx),'b-o');
        hBv = plot(w(lineIdx),Bvplot(lineIdx),'r-','linewidth',2);
        plotHandles = [hTotal hB hBv];
        plotNames = {'B+B_v (seakeeping)', 'B (potential)', ...
            'B_v(\omega) (viscous)'};

        if any(infIdx)
            hInf = plot(w(infIdx),Btotalplot(infIdx),'kx','linewidth',2);
            plot(w(infIdx),Bplot(infIdx),'bx','linewidth',2)
            plot(w(infIdx),Bvplot(infIdx),'rx','linewidth',2)
            plotHandles(end+1) = hInf;
            plotNames{end+1} = '\omega=\infty (stored at 10 rad/s)';
        end

        if hasPowerBeq
            Beqplot = B_eq(i,j)*ones(Nfreq,1);
            hBeq = plot(w(lineIdx),Beqplot(lineIdx),'c-.','linewidth',2);
            plotHandles(end+1) = hBeq;
            plotNames{end+1} = 'B_{eq} (power-based)';
            if any(infIdx)
                plot(w(infIdx),Beqplot(infIdx),'cx','linewidth',2)
            end
        end

        if hasPowerBeq && hasPowerBv
            D_eq = B_eq(i,j) + Bv_eq(i,j);
            Deqplot = D_eq*ones(Nfreq,1);
            hDeq = plot(w(lineIdx),Deqplot(lineIdx),'m-.','linewidth',2);
            plotHandles(end+1) = hDeq;
            plotNames{end+1} = 'B_{eq}+B_v (power-based)';
            if any(infIdx)
                plot(w(infIdx),Deqplot(infIdx),'mx','linewidth',2)
            end
        end

        hold off
        grid
        legend(plotHandles,plotNames,'Location','best','FontSize',8)
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');
        xlabel('frequency (rad/s)')
        title(Hw);
        k = k +1;        
    end
end
