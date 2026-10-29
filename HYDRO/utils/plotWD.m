function plotWD(vessel,x_axis,velno)
% plotWD plots signed mean wave-drift coefficients in surge, sway and yaw.
% Uses the current vessel.drift.w/coeff structure.
% Coefficients are loads per squared incident-wave amplitude. Moments retain
% the reference origin of the supplied dataset.
%
%    plotWD(vessel,x_axis,velno)
%
%    plotWD(vessel,'s',1)      - Period in seconds
%    plotWD(vessel,'hz',1)     - 1/s
%    plotWD(vessel,'rads',1)   - rad/s
%
% Inputs:
%    vessel: MSS vessel structure
%    x_axis {'rads','s','hz'}
%    velno  (optionally): speed number
%
% Author:    Thor I. Fossen
% Date:      2007-06-27
% Revisions:

if nargin < 2
    x_axis = 'rads';
end
if nargin < 3
    velno = 1;
elseif nargin > 3
    disp('Error: wrong number of input arguments')
    return
end

if isfield(vessel,'drift') && isstruct(vessel.drift) && ...
        isfield(vessel.drift,'w') && isfield(vessel.drift,'coeff')
    w = vessel.drift.w;
    amp = vessel.drift.coeff;
else
    error('The vessel has no vessel.drift mean wave-drift coefficients.');
end
validateattributes(w, {'numeric'}, {'vector','real','finite','positive'});
validateattributes(velno, {'numeric'}, {'scalar','integer','positive'});
if ~iscell(amp) || numel(amp) ~= 3
    error('vessel.drift.coeff must contain three cells for surge, sway and yaw.');
end
for j = 1:3
    values = amp{j};
    if ~isnumeric(values) || ~isreal(values) || any(~isfinite(values(:))) || ...
            size(values,1) ~= numel(w) || size(values,2) ~= 36 || ...
            velno > size(values,3)
        error('Invalid vessel.drift.coeff dimensions, values or speed index.');
    end
end
if ~any(strcmp(x_axis, {'rads','s','hz'}))
    error('x_axis must be rads, s or hz.');
end
txt   = 'Signed mean wave-drift coefficient:';
figno = 100;
ytxt11 = 'kN/m^2';
ytxt12 = 'kNm/m^2';
S1 = 1/1000;     % amp (N)  1-3 scaling
S2 = 1/1000;     % amp (Nm) 4-6 scaling

dim = size(amp{1});
if dim(2) ~= 36
    disp('Error: can only plot RAO data for 0-360 deg with 10 deg interval')
    return
end

% Surge
figure(figno)
arg = S1*amp{1}(:,:,velno);

if strcmp(x_axis,'rads')
   xarg = w;
   xtxt = 'Wave frequency (rad/s)';
elseif strcmp(x_axis,'s')
    xarg = 2*pi./w;
    xtxt = 'Period (s)';
else
    xarg = w/(2*pi);
    xtxt = 'Wave frequency (Hz)';
end

subplot(311)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SURGE'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid

% Sway
arg = S1*amp{2}(:,:,velno);
subplot(312)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SWAY'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid

% Yaw
subplot(313)
arg = S2*amp{3}(:,:,velno);
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' YAW'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt12)
grid

set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'type','legend'),'FontSize',10)
set(findall(gcf,'type','line'),'linewidth',1.5)

