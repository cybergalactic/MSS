function vessel = veres2vessel(filename, plot_flag)
% VERES2VESSEL (MSS Hydro)
%
% vessel = veres2vessel(filename, disp_flag) reads data from the
% ShipX (Veres) output files *.re1, *.re2, *.re7, *.re8, and *.hyd and
% store the data in vesselname.mat using the MSS vessel struture. Examples:
%
% >> veres2vessel('input') 
% >> veres2vessel('input','1111') 
%
% where input are the Veres output data file name.
% -------------------------------------------------------------------------
% Inputs:
%   filename: *  (without extension) reads and processes the following
%                 ShipX (Veres) files:                
%                    *.re1  motion RAOs (vessel.motionRAO)
%                    *.re2  wave drift data (vessel.driftfrc)
%                    *.re7  added mass, damping, restoring forces
%                    *.re8  force RAOsu (vessel.forceRAO)
%
%    plot_flag (optionally): '1000' plot A and B matrices
%			                 '0100' plot force RAOs
%                            '0010' plot motion RAOs
%                            '0001' plot wave drift forces
%                            '0000' NO PLOT
%                            '1111' PLOT ALL
% -------------------------------------------------------------------------
% Outputs:
% vessel contains data in the following form (saved to vesselname.mat):
%
%    vessel.headings:                     headings
%    vessel.velocities:                   velocities
%    vessel.freqs:                        frequencies (A and B matrices)
%
%    vessel.A(6,6,freqno,velno):          added mass matrix
%    vessel.B(6,6,freqno,velno):          damping matrix
%    vessel.C(6,6,freqno,velno):          spring stiffness matrix
%    vessel.MRB(6,6):                     rigid-body mass matrix
%
%    vessel.driftfrc.
%       amp(freqno,headno,velno):         wave drift force amplitudes
%       w(1,freqno):                      circular wave frequencies
%
%    vessel.forceRAO.
%       amp{1:6}(freqno,headno,velno):    wave excitation force amplitudes
%       phase{1:6}(freqno,headno,velno):  wave excitation force phases
%       w(1,freqno):                      circular wave frequencies
%
%    vessel.motionRAO.
%       amp{1:6}(freqno-n,headno,velno):   wave motion RAO amplitudes
%       phase{1:6}(freqno-n,headno,velno): wave motion RAO phases
%       w(1,freqno-n):                     circular wave frequencies where n>=0
%                                          is chosen such that w_o > w_min
%                                          (typically w_min = 0.2 rad/s)
%
%    vessel.main.
%       name:  vessel name
%       m:     mass
%       rho    density of water
%       nabla: displacement volume
%       k44:   radius of inertia 
%       k55:   radius of inertia 
%       k66:   radius of inertia 
%       m:     mass
%       CG:    centre of gravity  [LCG 0 VCG] w.r.t. Lpp/2 and Keel Line
%       CB:    centre of buoyancy [LCB 0 VCB] w.r.t. Lpp/2 and Keel Line
%       Lpp:   length between the perpendiculars
%       Lwl:   length of water line
%       T:     draught (water line)    
%       B:     breadth of ship
%       C_B:  block coefficient
%       S:    area of wetted surface
%       GM_L: longitudinal metacentric height
%       GM_T: transverse metacentric height
%
% Author:    T. I. Fossen
% Date:      2005-05-10 
% Revisions: 2008-02-15 Minor bug fixes
%            2009-09-11 Using new viscous damping viscous.m
%            2013-07-09 Fixed: MARINTEK coordinate origin not in CO
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%%
if ~exist('plot_flag')
	plot_flag = '1000';
end

disp(' ');
disp('************************* MSS Hydro *******************************')
disp('vessel = veres2vessel(*) computes the MSS vessel structure from the')
disp('ShipX (VERES) output files *.reN (N=1,2,7,8). The results are stored')
disp('on the file *.mat.')
disp(' ');
disp('Author: Thor I. Fossen');
disp('*******************************************************************')

fileinput  = input(['Vessel name used for output file (default = ' filename '): '],'s');
disp('Reading ShipX (VERES) data files *.re1, *.re2, *.re7, *.re8, and *.hyd...')

if strcmp(fileinput,'')
    vesselfile = filename;
else
    vesselfile = fileinput;       % default value
end

%--------------------------------------------------------------------------
%% read Shipx (Veres) *.re* data files
%--------------------------------------------------------------------------
vessel = read_veres_ABC(strcat(filename,'.re7'),str2num(plot_flag(1)));
data1  = read_veres_TF(strcat(filename,'.re8'), 0);
data2  = read_veres_TF(strcat(filename,'.re1'), 0);
data3  = read_veres_WD(strcat(filename,'.re2'), 0);

vessel.main.name = vesselfile;
vessel.forceRAO  = data1.forceRAO;
vessel.motionRAO = data2.motionRAO;
vessel.driftfrc  = data3.driftfrc;

%--------------------------------------------------------------------------
%% read Shipx (Veres) hydrostatic data from *.hyd file
%--------------------------------------------------------------------------
fid1 = fopen(strcat(filename,'.hyd'));

abort = 0;
while feof(fid1) == 0
    txt = char(fgetl(fid1));

    if  strfind(txt,'KB')
        KB   = str2num(txt(50:length(txt)-1));
    end
    if  strfind(txt,'LCB')
        LCB  = str2num(txt(50:length(txt)-1));
        txt  = char(fgetl(fid1));
        txt  = char(fgetl(fid1));
        C_B  = str2num(txt(50:length(txt)));
    end
    if  strfind(txt,'GMl')
        GM_L   = str2num(txt(50:length(txt)-1));
        txt  = char(fgetl(fid1));
        GM_T   = str2num(txt(50:length(txt)-1));
    end
end

vessel.main.GM_L  = GM_L;
vessel.main.GM_T  = GM_T;
vessel.main.C_B   = C_B;
vessel.main.CB(1) = LCB-vessel.main.Lpp/2;
vessel.main.CB(2) = 0;
vessel.main.CB(3) = KB;

fclose(fid1);

% approximations
vessel.main.Lwl = vessel.main.Lpp; 
vessel.main.S   = vessel.main.B*(vessel.main.Lpp + 2*vessel.main.T);

%--------------------------------------------------------------------------
%% plots
%--------------------------------------------------------------------------
if strcmp(plot_flag(2),'1')
    plotTF(vessel,'force','rads',1)
end
if strcmp(plot_flag(3),'1')
    plotTF(vessel,'motion','rads',1)
end
if strcmp(plot_flag(4),'1')
    plotWD(vessel,'rads',1)
end
%--------------------------------------------------------------------------
%% save data to <vessel>.mat
%--------------------------------------------------------------------------
save(vesselfile,'vessel')
disp(['Simulink structure <vessel> has been saved to: ' vesselfile,'.mat']);
disp('-------------------------------------------------------------------')
