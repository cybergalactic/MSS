function plotABC(vessel,mtrx,velno)
% plotABC plots A, B, or C as 3x3 longitudinal block on top of
% 3x3 lateral block below, with an empty spacer row between them.
%
%   plotABC(vessel,mtrx)       plots all speeds
%   plotABC(vessel,mtrx,velno) plots one speed only
%
% Inputs:
%   vessel : MSS vessel structure
%   mtrx   : 'A' added mass
%            'B' potential damping
%            'C' restoring forces
%   velno  : optional speed index
%
% Longitudinal block indices: [1 3 5] x [1 3 5]
% Lateral block indices:      [2 4 6] x [2 4 6]
%
% Author: Thor I. Fossen
% Date:   2026-04-05

    % Select matrix
    switch upper(mtrx)
        case 'A'
            H = vessel.A;
            figBase = 100;
            letter = 'A';
        case 'B'
            H = vessel.B;
            figBase = 200;
            letter = 'B';
        case 'C'
            H = vessel.C;
            figBase = 300;
            letter = 'C';
        otherwise
            error('mtrx must be ''A'', ''B'', or ''C''.');
    end

    freqs      = vessel.freqs(:);
    velocities = vessel.velocities(:);

    nvel = length(velocities);

    % Speed selection
    if nargin < 3
        velList = 1:nvel;
    else
        validateattributes(velno, {'numeric'}, ...
            {'scalar','integer','>=',1,'<=',nvel}, mfilename, 'velno');
        velList = velno;
    end

    % Index sets
    longIdx = [1 3 5];
    latIdx  = [2 4 6];

    % Use tiledlayout if available, otherwise fall back to subplot
    useTiled = exist('tiledlayout','file') == 2;

    for kvel = 1:length(velList)
        iv = velList(kvel);

        figure(figBase + iv);
        clf;

        if useTiled
            tl = tiledlayout(7,3,'TileSpacing','loose','Padding','compact');
            title(tl, sprintf('%s matrix, U = %.3g m/s', letter, velocities(iv)));
        end

        % Labels above the two 3x3 blocks
        annotation('textbox',[0.22 0.955 0.56 0.03], ...
            'String',sprintf('Longitudinal: surge, heave, pitch (U = %.3g m/s)', velocities(iv)), ...
            'EdgeColor','none', ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',12);

        annotation('textbox',[0.32 0.47 0.42 0.03], ...
            'String',sprintf('Lateral: sway, roll, yaw (U = %.3g m/s)', velocities(iv)), ...
            'EdgeColor','none', ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',12);

        % --- Top: longitudinal 3x3 block ---
        plotBlock(H, freqs, iv, longIdx, letter, 1, useTiled);

        % --- Bottom: lateral 3x3 block ---
        plotBlock(H, freqs, iv, latIdx, letter, 5, useTiled);
    end
end

function plotBlock(H, freqs, velno, idx, letter, rowOffset, useTiled)
% Helper to plot one 3x3 block.
    % VERES uses omega = 10 rad/s to represent infinite frequency.
    infIdx = abs(freqs - 10) < 10*eps(10);
    lineIdx = ~infIdx;

    count = 0;
    for rr = 1:3
        for cc = 1:3
            count = count + 1;
            i = idx(rr);
            j = idx(cc);

            Hij = reshape(H(i,j,:,velno),[],1);

            if useTiled
                nexttile((rowOffset-1)*3 + count);
            else
                subplot(7,3,(rowOffset-1)*3 + count);
            end

            plot(freqs(lineIdx),Hij(lineIdx),'b-o','MarkerSize',3)
            if any(infIdx)
                hold on
                plot(freqs(infIdx),Hij(infIdx),'rx', ...
                    'MarkerSize',7,'LineWidth',1.5)
                hold off
            end
            grid on;

            title(sprintf('%s_{%d%d}', letter, i, j), 'Interpreter','tex');
            xlabel('Frequency (rad/s)');
        end
    end
end
