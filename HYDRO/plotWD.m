function plotWD(vessel,x_axis,velno)
% plotWD (MSS Hydro)
%
% >> plotWD(vessel,x_axis,velno)      plots the wave drift amplitudes
%                                     versus seconds,frequency and rad/s
% >> plotWD(vessel,'s',1)      period in seconds
% >> plotWD(vessel,'hz',1)     1/s
% >> plotWD(vessel,'rads',1)   rad/s
%
% Input: 
%    vessel: MSS vessel structure
%    x_axis {'rads','s','hz'}
%    velno  (optionally):  speed number
%
% Author:    Thor I. Fossen
% Date:      2007-06-27
% Revisions: 
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

velocities = vessel.velocities;

nvel  = length(velocities);

if nargin == 2
    velno = 1;
elseif nargin > 3
    disp('Error: wrong number of input arguments')
    return
end

w     = vessel.driftfrc.w;
amp   = vessel.driftfrc.amp;
txt   = 'WD amplitude:';
figno = 100;
ytxt11 = 'kN/A^2';
ytxt12 = 'kNm/A^2';
S1 = 1/1000;     % amp (N)  1-3 scaling
S2 = 1/1000;     % amp (Nm) 4-6 scaling

dim = size(amp{1});
if dim(2) ~= 36,
    disp('Error: can only plot RAO data for 0-360 deg with 10 deg interval')
    return
end


% surge
figure(figno)
arg = S1*amp{1}(:,:,velno);

if strcmp(x_axis,'rads')
   xarg = w;
   xtxt = 'wave frequency (rad/s)';
elseif strcmp(x_axis,'s')
    xarg = 2*pi./w;
    xtxt = 'period (s)';
else
    xarg = w/(2*pi);
    xtxt = 'wave frequency (Hz)';    
end

subplot(311)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SURGE'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid

% sway
arg = S1*amp{2}(:,:,velno);
subplot(312)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SWAY'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid

% yaw
subplot(313)
arg = S2*amp{3}(:,:,velno);
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' YAW'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt12)
grid

