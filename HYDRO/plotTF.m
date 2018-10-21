function plotTF(vessel,type,x_axis,velno)
% plotTF (MSS Hydro)
%
% >> plotTF(vessel,type,x_axis, velno) plots the motion or force RAO
%                                      transfer functions  versus frequency
% >> plotTF(vessel,'motion','rads',1)  rad/s
% >> plotTF(vessel,'force','s',1)      period in seconds
% >> plotTF(vessel,'force','hz',1)     1/s
%
% Input: 
%    vessel: MSS vessel structure
%    type   {'motion','force'}  motion or force RAO
%    x_axis {'rads','s','hz'}
%    velno  (optionally):  speed number
%
% Author:    Thor I. Fossen
% Date:      2005-11-26
% Revisions: 
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

velocities = vessel.velocities;

nvel  = length(velocities);

if nargin == 3
    velno = 1;
elseif nargin > 4
    disp('Error: wrong number of input arguments')
    return
end

if strcmp(type,'motion')
    w     = vessel.motionRAO.w;
    amp   = vessel.motionRAO.amp;
    phase = vessel.motionRAO.phase;
    txt   = 'Motion RAO amplitude and phase:';
    figno = 70;    
    ytxt11 = 'm';
    ytxt12 = 'deg';
    ytxt2  = 'deg';
    S1 = 1;           % amp (m)   1-3 scaling
    S2 = 180/pi;      % amp (rad) 4-6 scaling
else
    w     = vessel.forceRAO.w;
    amp   = vessel.forceRAO.amp;
    phase = vessel.forceRAO.phase;
    txt   = 'Force RAO amplitude and phase:';
    figno = 80;    
    ytxt11 = 'kN';
    ytxt12 = 'kNm';
    ytxt2  = 'deg';
    S1 = 1/1000;     % amp (N)  1-3 scaling
    S2 = 1/1000;     % amp (Nm) 4-6 scaling
end

dim = size(amp{1});
if dim(2) ~= 36
    disp('Error: can only plot RAO data for 0-360 deg with 10 deg interval')
    return
end


% surge
figure(figno)
arg = S1*amp{1}(:,:,velno);
phs = (180/pi)*unwrap(phase{1}(:,:,velno));

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

subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SURGE'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid

% sway
figure(figno+1)
arg = S1*amp{2}(:,:,velno);
phs = (180/pi)*unwrap(phase{2}(:,:,velno));

subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' SWAY'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid

% heave
figure(figno+2)
arg = S1*amp{3}(:,:,velno);
phs = (180/pi)*unwrap(phase{3}(:,:,velno));

subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' HEAVE'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt11)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid

%  roll
figure(figno+3)
arg = S2*amp{4}(:,:,velno);
phs = (180/pi)*unwrap(phase{4}(:,:,velno));

subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' ROLL'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt12)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid

% pitch
figure(figno+4)
arg = S2*amp{5}(:,:,velno);
phs = (180/pi)*unwrap(phase{5}(:,:,velno));

subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' PITCH'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt12)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid

figure(figno+5)
arg = S2*amp{6}(:,:,velno);
phs = (180/pi)*unwrap(phase{6}(:,:,velno));

% yaw
subplot(211)
plot(xarg,arg(:,1),xarg,arg(:,4),xarg,arg(:,7),xarg,arg(:,10),xarg,arg(:,13),xarg,arg(:,16),xarg,arg(:,19))
title(strcat(txt,' YAW'))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt12)
grid
subplot(212)
plot(xarg,phs(:,1),xarg,phs(:,4),xarg,phs(:,7),xarg,phs(:,10),xarg,phs(:,13),xarg,phs(:,16),xarg,phs(:,19))
legend('0 deg','30 deg','60 deg','90 deg','120 deg','150 deg','180 deg')
xlabel(xtxt)
ylabel(ytxt2)
grid
