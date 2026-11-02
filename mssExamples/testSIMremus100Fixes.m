function testSIMremus100Fixes()
% testSIMremus100Fixes() runs the first 600 s of SIMremus100.m with control
% method 3 (2-D ALOS, depth reference model) and checks two things.
%
% 1. Yaw acceleration. In method 3 the heading reference comes from
%    LOSobserver, which produces psi_d and r_d and no yaw acceleration, so
%    the a_d that integralSMCheading reads (a yaw acceleration, rad/s^2) must
%    stay 0. The depth reference model refModel(z_d,zdot_d,a_d,...) of the
%    same method used the same variable a_d for its vertical acceleration
%    (m/s^2), so after the first change of the depth reference the
%    heading controller received a vertical acceleration as feed-forward.
%
% 2. Propeller ramp. The ramp raises n_p from 1000 rpm by n_rate = 0.1 rpm per
%    step towards n_d = 1300 rpm. After it has arrived (150 s) n_p must stay
%    at n_d; the two-sided test "if n_p < n_d ... elseif n_p > n_d" with
%    floating-point sums makes it alternate between 1300.0 and 1300.1 rpm.
%
% The test does not change SIMremus100.m: it runs a copy in which the dialog
% controlMethod() is replaced by the setting "method 3, Euler angles", and it
% records the arguments passed to integralSMCheading and the propeller speed
% passed to remus100 through two thin wrapper functions placed first on the
% path. The run is stopped after 12000 steps (600 s).
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Reference:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%       Control, 3rd ed. Section 12.1.1, Eq. 12.10, the reference model of
%       refModel.m (as cited in that file); a_d is the acceleration of the
%       quantity the model filters, here the depth in m.
%
% Author:    Enio Krizman
% Date:      2 November 2026

nStop = 12000;                       % integralSMCheading calls (steps)
n_d = 1300;                          % propeller speed target in SIMremus100.m

tmp = tempname; mkdir(tmp);
cleanupAll = onCleanup(@() removeTemp(tmp));
set(0,'DefaultFigureVisible','off');

realRemus = @remus100;
realSMC   = @integralSMCheading;
setappdata(0,'mq_real_remus100',realRemus);
setappdata(0,'mq_real_smc',realSMC);
setappdata(0,'mq_n_log',zeros(1,4*nStop+16));
setappdata(0,'mq_n_cnt',0);
setappdata(0,'mq_a_log',zeros(1,nStop+4));
setappdata(0,'mq_a_cnt',0);
setappdata(0,'mq_nStop',nStop);

src = fileread(which('SIMremus100'));
src = strrep(src,'[ControlFlag, KinematicsFlag] = controlMethod();', ...
    'ControlFlag = 3; KinematicsFlag = 1;');
src = strrep(src,'function SIMremus100()','function SIMremus100_mq()');
fid = fopen(fullfile(tmp,'SIMremus100_mq.m'),'w'); fwrite(fid,src); fclose(fid);

fid = fopen(fullfile(tmp,'remus100.m'),'w');
fprintf(fid,['function varargout = remus100(varargin)\n' ...
    'real = getappdata(0,''mq_real_remus100'');\n' ...
    'if numel(varargin) >= 2\n' ...
    '  k = getappdata(0,''mq_n_cnt'') + 1; L = getappdata(0,''mq_n_log'');\n' ...
    '  L(k) = varargin{2}(3); setappdata(0,''mq_n_log'',L); setappdata(0,''mq_n_cnt'',k);\n' ...
    'end\n' ...
    '[varargout{1:max(nargout,1)}] = real(varargin{:});\n' ...
    'end\n']);
fclose(fid);

fid = fopen(fullfile(tmp,'integralSMCheading.m'),'w');
fprintf(fid,['function varargout = integralSMCheading(varargin)\n' ...
    'real = getappdata(0,''mq_real_smc'');\n' ...
    'k = getappdata(0,''mq_a_cnt'') + 1; L = getappdata(0,''mq_a_log'');\n' ...
    'L(k) = varargin{5}; setappdata(0,''mq_a_log'',L); setappdata(0,''mq_a_cnt'',k);\n' ...
    'if k > getappdata(0,''mq_nStop''), error(''MQ:stop'',''stop''); end\n' ...
    '[varargout{1:max(nargout,1)}] = real(varargin{:});\n' ...
    'end\n']);
fclose(fid);

addpath(tmp);
stopped = false;
try
    SIMremus100_mq();
catch ME
    if strcmp(ME.identifier,'MQ:stop'), stopped = true; else, rethrow(ME); end
end
if ~stopped
    error('testSIMremus100Fixes: the simulation did not reach the stop step');
end

aLog = getappdata(0,'mq_a_log'); aLog = aLog(1:getappdata(0,'mq_a_cnt')-1);
nLog = getappdata(0,'mq_n_log'); nLog = nLog(1:getappdata(0,'mq_n_cnt'));
nfail = 0;

aMax = max(abs(aLog));
if aMax > 0
    fprintf('yaw acceleration a_d read by integralSMCheading   FAILED: max |a_d| = %.3e (must be 0)\n', aMax);
    nfail = nfail + 1;
else
    fprintf('yaw acceleration a_d read by integralSMCheading   PASSED: max |a_d| = 0\n');
end

late = nLog(end-4*2000+1:end);       % last 2000 steps, ramp finished at 3000
spread = max(late) - min(late);
miss = max(abs(late - n_d));
if spread > 1e-6 || miss > 1e-6
    fprintf('propeller speed after the ramp                    FAILED: spread %.3g rpm, max |n_p - n_d| %.3g rpm\n', spread, miss);
    nfail = nfail + 1;
else
    fprintf('propeller speed after the ramp                    PASSED: spread %.3g rpm\n', spread);
end

if nfail > 0
    error('testSIMremus100Fixes FAILED (%d check(s) wrong)', nfail);
end
disp('testSIMremus100Fixes PASSED');
end

function removeTemp(tmp)
warnState = warning('off','all');
rmpath(tmp); rmdir(tmp,'s');
warning(warnState);
cleanAppdata();
end

function cleanAppdata()
for f = {'mq_real_remus100','mq_real_smc','mq_n_log','mq_n_cnt','mq_a_log','mq_a_cnt','mq_nStop'}
    if isappdata(0,f{1}), rmappdata(0,f{1}); end
end
end
