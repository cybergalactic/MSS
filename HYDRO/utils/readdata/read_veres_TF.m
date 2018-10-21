function vessel = read_veres_TF(filename, disp_flag)
% vessel = read_veres_TF(filename, disp_flag)
%
% Read data from Veres output file *.re1 or *.re8
%
% Inputs:
%
% filename: 	Text string with name of VERES output file
% disp_flag: 	1 for displaying results to screen
%				2 for also displaying TF amplitude plots
%				3 for adding TF phase plots
%				0 for none of the above
%
% Outputs:
%
% TF_data contains transfer function data extracted from the file in form
% of amplitude amplification and phase, structured as:
%
% 	vessel.forceRAO.amp{dofno}(freqno,headnno,velno)   
% 	vessel.forceRAO.phase{dofno}(freqno,headnno,velno) 
%
% 	vessel.motionRAO.amp{dofno}(freqno,headnno,velno)  
% 	vessel.motionRAO.phase{dofno}(freqno,headnno,velno) 
% 	vessel.motionRAO.w(1,freqno)
%
% 	vessel.headings   	
% 	vessel.velocities	
% 
% 'vel', 'dir' and 'freq' contain velocities in m/s, directions in radians
% and frequencies in rad/s respectively.
%
% 'Amp' and 'Phase' are 3D matrices with dimensions (freq, phase, vel)
%
% Author:    ?. N. Smogeli and Thor I. Fossen
% Date:      2004-04-12 Beta version
% Revisions: 2007-12-10 T. I. Fossen - R1.0
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

if ~exist('disp_flag')
	disp_flag = 3;
end

%% ------------------------------------------------------------------------
% Read data from file
%--------------------------------------------------------------------------

% Open file
fid = fopen(filename);

if fid <= 0
	disp('Invalid filename')
end

% Read the first 6 lines, which contain written run info
for m = 1:6	
	header{m} = fgetl(fid);	
end

% Read the next 4 lines, which contain numeric run info
for m = 7:10
	header{m} = fgets(fid);
end

% Extract vessel info
temp = str2num(header{7});
rho = temp(1);
g   = temp(2);
temp = str2num(header{8});
LPP = temp(1);
B = temp(2);
T = temp(3);
temp = str2num(header{9});
LCG = temp(1);
VCG = temp(2);

vessel.main.rho = rho;
vessel.main.Lpp = LPP;
vessel.main.B   = B;
vessel.main.T   = T;
vessel.main.CG  = [-LCG 0 VCG];   % x-forward

% Extract run info from line 10
temp = str2num(header{10});

nvel = temp(1);
nhead = temp(2);
nfreq = temp(3);
ndof = temp(4);

% Loop over velocities
for velno = 1:nvel

    k = 0; % Response index
	
    vels{velno} = str2num(fgets(fid));
	
	% Loop over headings
	for headno = 1:nhead
		
		headings(headno) = str2num(fgets(fid));
		
		% Loop over frequencies
		for freqno = 1:nfreq
			
			freqs(freqno) = str2num(fgets(fid));
			
			% Update response index
			k = k + 1;
			
			% Loop over DOF
			for dofno = 1:ndof
				
				% Read [DOF Real Imag], store as amplitude and phase
				temp = str2num(fgets(fid));
				resp{velno}.amp(k,temp(1)) = abs(temp(2) + temp(3)*i); 
				resp{velno}.phase(k,temp(1)) = angle(temp(2) + temp(3)*i); 
				
			end % Loop over DOF
			
		end % Loop over frequencies
		
	end % Loop over headings
	
end % Loop over velocities

% Extract velocity info
for m = 1:nvel
	vel(m) = vels{m}(1);
	sink(m) = vels{m}(2);
	trim(m) = vels{m}(3);
	xmtn(m) = vels{m}(4);
	zmtn(m) = vels{m}(5);	
end

% check to see if raw data for 10:10:180 deg is available
for ii=1:19
    if headings(ii) ~= 10*(ii-1);
        disp('Error: MSS Hydro requires RAO data for headings 0:10:180 deg');
        return
    end
end
%% ------------------------------------------------------------------------
% Scale, transform, and store data in MSS vessel structure
% (x-forward, y-starboard, z-upwards, 0 deg beam seas)
%
% A_new   = T*A*T,   Veres to MSS (Fossen 2002): T = diag([-1 1 -1 -1 1 -1])
% tau_new = T*tau    Veres (0 deg head seas), MSS (0 deg beam seas)
%--------------------------------------------------------------------------
nhead_tot = (nhead - 2)*2 + 2;

% preallocate memory for maximum 10 velocties
ext = filename(length(filename)-2:length(filename));
if ext == 're8',
    for dofno = 1:6
        vessel.forceRAO.amp{dofno}   = zeros(nfreq,nhead_tot,10);
        vessel.forceRAO.phase{dofno} = zeros(nfreq,nhead_tot,10);
    end
elseif ext == 're1',
    for dofno = 1:6
        vessel.motionRAO.amp{dofno}   = zeros(nfreq,nhead_tot,10);
        vessel.motionRAO.phase{dofno} = zeros(nfreq,nhead_tot,10);
    end
end

for dofno = 1:ndof

	for velno = 1:nvel

		% Unwrap the phase
		resp{velno}.phase(:,dofno) = unwrap(resp{velno}.phase(:,dofno));
        
        % Get on matrix form
		amp_mat = reshape(resp{velno}.amp(:,dofno),nfreq,nhead);
		phs_mat = reshape(resp{velno}.phase(:,dofno),nfreq,nhead);
        
        % Change sequence to get the headings correct (in Veres the
        % headings are defined relative to the bow while the MSS standard is
        % relative to the stern with x-axis forward, i.e 180 deg difference).
        for head = 1:nhead         
            amp(:,head) = amp_mat(:,nhead + 1 - head);
            phs(:,head) = phs_mat(:,nhead + 1 - head);
        end
		
        % Correct phases in surge, heave, roll and yaw (1, 3, 4, 6)
        % corresponding to the matrix: T = diag([-1 1 -1 -1 1 -1])
        % Amplitudes are kept positive
        if dofno == 1 | dofno == 3 | dofno == 4 | dofno == 6
            phs = phs + pi;           % rad
        end
        
        %add data for 180-360 deg (symmetry)
        for ii = 2:(nhead-1)
            jj =  nhead_tot + 2 - ii;
            headings(jj) = -headings(ii) + 360;

            if dofno == 1 || dofno == 3 || dofno == 5
                amp(:,jj) = amp(:,ii);
                phs(:,jj) = phs(:,ii);
            else
                amp(:,jj) = amp(:,ii);
                phs(:,jj) = rad2pipi(phs(:,ii) - pi);
            end
        end
                  
        % fill in RAO data                   
        if ext == 're8',
            vessel.forceRAO.amp{dofno}(:,:,velno)   = amp;
    		vessel.forceRAO.phase{dofno}(:,:,velno) = phs;
     		vessel.forceRAO.w                       = freqs;
        elseif ext == 're1',          
            vessel.motionRAO.amp{dofno}(:,:,velno)   = amp;
    		vessel.motionRAO.phase{dofno}(:,:,velno) = phs;
     		vessel.motionRAO.w                       = freqs;            
        end
        
        vessel.headings   	= headings*pi/180;
		vessel.velocities	= vel;
		
	end
	
end

%% ------------------------------------------------------------------------
% Display run info on screen
%--------------------------------------------------------------------------
if disp_flag > 0
	
	disp(' ')
	disp('-------------------------------------------')
	disp('       VERES TRANSFER FUNCTION DATA')
	disp('-------------------------------------------')
	disp('File header: ')
	disp(header{2})
	disp(header{3})
	disp(header{4})
	disp(header{5})
	disp(header{6})
	disp('-------------------------------------------')
	disp('Run info:')
	disp(' Main particulars:')
	disp(sprintf(' LPP     : %0.2f m',LPP))
	disp(sprintf(' Breadth : %0.2f m',B))
	disp(sprintf(' Draft   : %0.2f m',T))
	disp(sprintf(' LCG     : %0.2f m (relative to LPP/2)',LCG))
	disp(sprintf(' VCG     : %0.2f m (relative to baseline)',VCG))
	disp(' ')
	disp(sprintf(' Number of velocities  : %d',nvel))
	disp(sprintf(' Number of headings    : %d',nhead_tot))
	disp(sprintf(' Number of frequencies : %d',nfreq))
	disp(sprintf(' Number of DOF         : %d',ndof))
	disp('-------------------------------------------')
	disp('Details:')
    txt = sprintf('%d ',headings);
    disp(['Headings (deg)     : ' ,txt]);
    txt = sprintf('%1.2f ',freqs);    
    disp(['Frequencies (rad/s): ',txt]);
    disp(' ')
	for m = 1:nvel
		disp(sprintf(' Velocity %d  : %0.2f m/s',m,vel(m)))
		disp(sprintf('    Sinkage  : %0.2f m',sink(m)))
		disp(sprintf('    Trim     : %0.2f deg',trim(m)))
		disp(sprintf('    Motion defined in:'))
		disp(sprintf('    X-Coord  : %0.2f m (rel. LPP/2)',xmtn(m)))
		disp(sprintf('    Y-Coord  : 0.00 m (rel. centerline)'))
		disp(sprintf('    Z-Coord  : %0.2f m (rel. baseline)',zmtn(m)))
		disp(' ')
	end
	disp('-------------------------------------------')
	
end

%% ------------------------------------------------------------------------
% Plot transfer function
%--------------------------------------------------------------------------
if disp_flag > 1

    if ext == 're8',
        RAO = vessel.forceRAO;
    elseif ext == 're1',
        RAO = vessel.motionRAO;
    end
    
	for velno = 1:nvel
        if  any(any(RAO.amp{2}(:,:,velno)))
            for dofno = 1:ndof
                amp = RAO.amp{dofno}(:,:,velno);
                phs = RAO.phase{dofno}(:,:,velno);
                figure
                mesh(headings, freqs, amp)
                xlabel('Directions [deg]')
                ylabel('Frequencies [rad/s]')
                title(sprintf('RAO amplitude [-] for DOF %d and V = %0.1f m/s',dofno,vels{velno}(1)))

                if disp_flag > 2
                    figure
                    mesh(headings, freqs, phs)
                    xlabel('Directions [deg]')
                    ylabel('Frequencies [rad/s]')
                    title(sprintf('RAO phase [rad] for DOF %d and V = %0.1f m/s',dofno,vels{velno}(1)))
                end
            end
        end
	end
	
end