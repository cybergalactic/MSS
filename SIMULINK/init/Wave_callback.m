function void = Wave_callback
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% void = Wave_callback
%
% Project:	"Marine System Simulator: MSS"
%
% Abstract:	Enable / disable mask parameters for different spectra in Simulink Waves block
%
% Inputs:	-
%
% Outputs:	-
%
% Calls:	-
%
% Author:	Oyvind Smogeli, March 2005
%
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
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

% Get the current value of parameter spectrum_type
choice = get_param(gcb,'spectrum_type');

% Find the possible choices for spectrum_type
enum = get_param(gcb,'Objectparameters');
enum = enum.spectrum_type.Enum;

% Convert to number
for q = 1:length(enum)
	if strcmp(choice, enum{q})
		type = q;
	end	
end

on = 'on';
off = 'off';

% Load current settings
onoff = get_param(gcb,'MaskVisibilities');

% Cell array of dialog parameters
diag_params = fieldnames(get_param(gcb,'DialogParameters'));

% Number of parameters
nparams = length(diag_params);

% Initally set all parameters disabled
for i = 1:nparams
	onoff{i} = off;
end
	
% counter 
m = 1;

% Set the common fields enabled
on_id{m} = 'spectrum_type'; m=m+1;

% Set the desired parameters for the spectrum type active
% 1 = ITTC, 2 = JONSWAP, 3 = Torsethaugen Doubly Peaked, 4 = predefined waves
switch type
	
	case 1
		
		% Enabled parameters for ITTC
		on_id{m} = 'hs'; m=m+1;
		on_id{m} = 'omega_peak'; m=m+1;
		on_id{m} = 'psi_mean'; m=m+1;
		on_id{m} = 'spread'; m=m+1;
		on_id{m} = 'depth'; m=m+1;
		on_id{m} = 'nfreq'; m=m+1;
		on_id{m} = 'ndir'; m=m+1;
		on_id{m} = 'energylim'; m=m+1;
		on_id{m} = 'freq_cutoff'; m=m+1;
		on_id{m} = 'dir_cutoff'; m=m+1;
		on_id{m} = 'rand_freq'; m=m+1;
		on_id{m} = 'rand_dir'; m=m+1;
		on_id{m} = 'rand_seed'; m=m+1;
		on_id{m} = 'plot_spectrum'; m=m+1;
		on_id{m} = 'plot_realization'; m=m+1;
		on_id{m} = 'disp_flag'; m=m+1;
		
	case 2
		
		% Enabled parameters for JONSWAP
		on_id{m} = 'hs'; m=m+1;
		on_id{m} = 'omega_peak'; m=m+1;
		on_id{m} = 'psi_mean'; m=m+1;
		on_id{m} = 'spread'; m=m+1;
		on_id{m} = 'depth'; m=m+1;
		on_id{m} = 'nfreq'; m=m+1;
		on_id{m} = 'ndir'; m=m+1;
		on_id{m} = 'energylim'; m=m+1;
		on_id{m} = 'freq_cutoff'; m=m+1;
		on_id{m} = 'dir_cutoff'; m=m+1;
		on_id{m} = 'rand_freq'; m=m+1;
		on_id{m} = 'rand_dir'; m=m+1;
		on_id{m} = 'rand_seed'; m=m+1;
		on_id{m} = 'plot_spectrum'; m=m+1;
		on_id{m} = 'plot_realization'; m=m+1;
		on_id{m} = 'disp_flag'; m=m+1;
		on_id{m} = 'gamma'; m=m+1;
	
	case 3
		
		% Enabled parameters for Torsethaugen
		on_id{m} = 'hs'; m=m+1;
		on_id{m} = 'omega_peak'; m=m+1;
		on_id{m} = 'psi_mean'; m=m+1;
		on_id{m} = 'spread'; m=m+1;
		on_id{m} = 'depth'; m=m+1;
		on_id{m} = 'nfreq'; m=m+1;
		on_id{m} = 'ndir'; m=m+1;
		on_id{m} = 'energylim'; m=m+1;
		on_id{m} = 'freq_cutoff'; m=m+1;
		on_id{m} = 'dir_cutoff'; m=m+1;
		on_id{m} = 'rand_freq'; m=m+1;
		on_id{m} = 'rand_dir'; m=m+1;
		on_id{m} = 'rand_seed'; m=m+1;
		on_id{m} = 'plot_spectrum'; m=m+1;
		on_id{m} = 'plot_realization'; m=m+1;
		on_id{m} = 'disp_flag'; m=m+1;

	case 4
		
		% Enabled parameters for predefined waves
		on_id{m} = 'Zeta_a'; m=m+1;
		on_id{m} = 'Omega'; m=m+1;
		on_id{m} = 'Phase'; m=m+1;
		on_id{m} = 'Wavenum'; m=m+1;
		on_id{m} = 'Psi'; m=m+1;
		
	otherwise
		
		disp('Error in Wave_callback: Undefined spectrum type!')
		
end
	
% Find the index of the parameters to be turned on
for i = 1:(m-1)

	on_index(i) = find(strcmp(on_id(i),diag_params) == 1);

end

% Set active values
for i = on_index
	
	onoff{i} = on;
	
end

% Modify the dialog 
set_param(gcb,'MaskVisibilities',onoff);

