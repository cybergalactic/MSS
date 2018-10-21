function vessel = read_veres_WD(filename, disp_flag)
% vessel = read_veres_WD(filename, disp_flag)
%
% Read and process data from Veres output file *.re2
%
% Inputs:
%   filename: 	Text string *.re2 with name of VERES output file
%
%   disp_flag: 	1 for displaying results to screen
%				2 for also displaying WD amplitude plots
%				0 for none of the above
%
% Outputs (MSS vessel structure):
%
%   vessel.driftfrc.
%       amp{dof}(freqno,headno,velno):    wave drift force amplitudes for
%                                         dof = 1,2,6
%       w(1,freqno):                      circular wave frequencies
%
%		vessel.headings:                  wave directions
%		vessel.velocities:                vessel speeds
%
% For a given wave, the wavedrift force F is given as:
%
%     F = vessel.driftfrc.amp{dof}*zeta_a^2,   dof = 1,2,6
%
% where amp{dof} is the amplitude magnification found from the wavedrift table
% and zeta_a is the wave amplitude.
%
% Authors:   ?. N. Smogeli and T. I. Fossen
% Date:      2005-05-11
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

% Read the next 3 lines, which contain numeric run info
for m = 7:9
	header{m} = fgets(fid);
end

% Extract vessel info
temp = str2num(header{7});
rho = temp(1);
g = temp(2);

temp = str2num(header{8});
LPP = temp(1);
B = temp(2);
T = temp(3);

% Extract run info from line 9
temp = str2num(header{9});
nvel = temp(1);
nhead = temp(2);
nfreq = temp(3);

% Loop over velocities
for velno = 1:nvel
    
    k = 0;     % Response index
	
	% Loop over headings
	for headno = 1:nhead

		temp = str2num(fgets(fid));
		vels(velno) = temp(1);
		headings(headno) = temp(2);
		
		% Loop over frequencies
		for freqno = 1:nfreq
			
			temp = str2num(fgets(fid));
			
			freqs(freqno) = temp(1);
						
			% Update response index
			k = k + 1;
			
			% store amplitude coef for surge, sway, yaw
			resp{velno}.amp(k,1) = temp(2)*rho*g*B^2/LPP;
			resp{velno}.amp(k,2) = temp(3)*rho*g*B^2/LPP;
			resp{velno}.amp(k,3) = temp(4)*rho*g*B^2;
			
		end % Loop over frequencies
		
	end % Loop over headings
	
end % Loop over velocities

%close file
fclose(fid);

% check to see if raw data for 10:10:180 deg is available
for ii=1:19
    if headings(ii) ~= 10*(ii-1);
        disp('Error: MSS Hydro requires RAO data for headings 0:10:180 deg');
        return
    end
end

%% ------------------------------------------------------------------------
% Transform the data to the h-frame
%--------------------------------------------------------------------------
nhead_tot = (nhead - 2)*2 + 2;

% preallocate memory for maximum 10 velocties
for dofno = 1:6
    vessel.driftfrc.amp{dofno} = zeros(nfreq,nhead_tot,10);
end

for dofno = 1:3

	for velno = 1:nvel
			
		amp_mat = reshape(resp{velno}.amp(:,dofno),nfreq,nhead);
		
		% Change sequence to get the headings correct (in Veres the
		% headings are defined relative to the bow, in MSS relative
		% to the stern with x-axis forward).
        for head = 1:nhead
            amp(:,head) = amp_mat(:,nhead + 1 - head);
        end
		
        % Change signs of amplitudes in surge and yaw (1, 3) - MSS axes
        if dofno == 1 | dofno == 3
            amp = -amp;    
        end

        %add data for 180-360 deg (symmetry)
        for ii = 2:(nhead-1)
            jj =  nhead_tot + 2 - ii;
            headings(jj) = -headings(ii) + 360;

            if dofno == 1 || dofno == 3 || dofno == 5
                amp(:,jj) = amp(:,ii);
            else
                amp(:,jj) = -amp(:,ii);
            end
        end
               
        % fill in RAO data
		vessel.driftfrc.amp{dofno}(:,:,velno) = amp;
		vessel.driftfrc.w	= freqs;
        
		vessel.headings	    = headings*pi/180;
		vessel.velocities	= vels;
				
	end
	
end

%% ------------------------------------------------------------------------
% Display run info on screen
%--------------------------------------------------------------------------
if disp_flag > 0
	
	disp(' ')
	disp('-------------------------------------------')
	disp('       VERES WAVEDRIFT DATA')
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
	disp(' ')
	disp(sprintf(' Number of velocities  : %d',nvel))
	disp(sprintf(' Number of headings    : %d',nhead_tot))
	disp(sprintf(' Number of frequencies : %d',nfreq))
	disp('-------------------------------------------')
	disp('Details:')
    txt = sprintf('%d ',headings);
    disp(['Headings (deg)     : ' ,txt]);
    txt = sprintf('%1.2f ',freqs);    
    disp(['Frequencies (rad/s): ',txt]);
	disp(['Velocities (m/S)   : ', num2str(vels)])
	disp('-------------------------------------------')
	
end

%% ------------------------------------------------------------------------
% Plot data
%--------------------------------------------------------------------------
if disp_flag > 1		
	for velno = 1:nvel		
		for dofno = 1:3	
			amp = vessel.driftfrc.amp{dofno}(:,:,velno);
			figure
			mesh(headings, freqs, amp)
			xlabel('Directions [deg]')
			ylabel('Frequencies [rad/s]')
			title(sprintf('WD amplitude [-] for DOF %d and V = %0.1f m/s',dofno,vels(velno)))	
		end
    end	
end