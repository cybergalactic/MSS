function plotRAO(vessel,headings,blackWhite)
% Function for plotting motion and force RAOs. For a list of available MSS
% vessel names with MAT-files containing the structure 'vessel', run
%
%   mssCraft
%
% Example calls:
%   load 'supply' vessel
%   plotRAO(vessel)                     % Default: [0 30 60 90], colored lines
%   plotRAO(vessel,[0 30 60 90])        % Selected headings, colored lines
%   plotRAO(vessel,[0 30 60 90],false)  % Same as the default
%   plotRAO(vessel,0:30:180,true)       % Seven headings, black-and-white lines
%   plotRAO(vessel,[],true)             % Default headings, black-and-white lines
%
% headings is a vector in degrees, matched to vessel.headings (radians).
% Requested headings must be present in the dataset; no interpolation is used.
% Omitted or empty headings defaults to [0 30 60 90]. blackWhite is an
% optional logical or 0/1 flag and defaults to false (colored lines).
% Figures 1 and 2 show motion RAOs; figures 3 and 4 show force RAOs.
% Each figure has three rows (one per DOF), with amplitude on the left
% and phase on the right. The first speed case is used.
%
% Author:    Thor I. Fossen
% Date:      2026-10-27
% Revisions:

if nargin < 2 || isempty(headings)
    headings = [0 30 60 90];
end
if nargin < 3
    blackWhite = false;
end
validateattributes(headings, {'numeric'}, ...
    {'vector','real','finite','nonempty'}, mfilename, 'headings');
columns = zeros(1,numel(headings));
labels = cell(1,numel(headings));
for i = 1:numel(headings)
    distance = abs(mod(vessel.headings(:)*180/pi - headings(i) + 180,360) - 180);
    [difference, columns(i)] = min(distance);
    if isempty(difference) || difference > 1e-6
        error('plotRAO:HeadingNotAvailable', ...
            'Requested heading %g deg is not available in vessel.headings.', headings(i));
    end
    labels{i} = sprintf('%g deg',headings(i));
end
validateattributes(blackWhite, {'logical','numeric'}, ...
    {'scalar','real','finite','binary'}, mfilename, 'blackWhite');

pre  = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
post = {' (m)',' (m)',' (m)',' (deg)',' (deg)',' (deg)'};
txtMa = strcat(pre,' Motion RAO - Amplitude',post); 
txtMp = strcat(pre,' Motion RAO - Phase (deg)'); 
post = {' (N)',' (N)',' (N)',' (Nm)',' (Nm)',' (Nm)'};
txtFa = strcat(pre,' Force RAO - Amplitude',post); 
txtFp = strcat(pre,' Force RAO - Phase (deg)'); 

for figureNo = 1:4
    figure(figureNo);
    clf;
    if figureNo <= 2
        rao = vessel.motionRAO;
        amplitudeTitles = txtMa;
        phaseTitles = txtMp;
    else
        rao = vessel.forceRAO;
        amplitudeTitles = txtFa;
        phaseTitles = txtFp;
    end

    firstDOF = 1 + 3*mod(figureNo - 1,2);
    for row = 1:3
        DOF = firstDOF + row - 1;
        subplot(3,2,2*row - 1);
        plotRAOamp(rao.w,rao.amp,DOF,figureNo <= 2,blackWhite,columns,labels);
        title(amplitudeTitles{DOF});
        subplot(3,2,2*row);
        plotRAOphs(rao.w,rao.phase,DOF,blackWhite,columns,labels);
        title(phaseTitles{DOF});
    end

    set(findall(gcf,'type','line'),'linewidth',1);
    set(findall(gcf,'type','text'),'FontSize',12);
    set(findall(gcf,'type','legend'),'FontSize',10);
end

end

%% FUNCTIONS
function plotRAOamp(w,amp,DOF,isMotion,blackWhite,columns,labels)

    velno = 1;   
    arg = amp{DOF}(:,:,velno);

    % Scale angular motion amplitudes to degrees; force moments stay in Nm.
    S = 1; if isMotion && DOF > 3
        S = 180/pi; 
    end  
    plotHeadings(w,S*arg(:,columns),blackWhite,labels);
end

function plotRAOphs(w,phase,DOF,blackWhite,columns,labels)

    velno = 1;    
    phs = (180/pi)*unwrap(phase{DOF}(:,:,velno));  
    plotHeadings(w,phs(:,columns),blackWhite,labels);
end

function plotHeadings(w,values,blackWhite,labels)
    lines = plot(w,values);
    markers = {'*','o','<','x','s','d','^','v','>','+','p','h'};
    styles = {'-','--',':','-.'};
    for i = 1:numel(lines)
        set(lines(i),'Marker',markers{1 + mod(i-1,numel(markers))}, ...
            'LineStyle',styles{1 + mod(floor((i-1)/numel(markers)),numel(styles))});
    end
    if blackWhite
        set(lines,'Color','k');
    end
    legend(lines,labels);
    xlabel('Frequency (rad/s)');
    grid on
end
