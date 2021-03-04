function void = MSS_mask_page(npar)
% Set the visibility of simulink parameters dialog with npar parameter per page
% The first parameter of the dialog must be called page and be a popup menu
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

% Get the current value of parameter page
page = str2num(get_param(gcb,'page'));

% Number of dialog parameters
n = length(fieldnames(get_param(gcb,'DialogParameters')));

on = 'on';
off = 'off';

% Load current settings
onoff = get_param(gcb,'MaskVisibilities');

% Initialize values
onoff{1} = on;

for i = 2:n

	onoff{i} = off;
end

% Set active value

for i = ((page-1)*npar + 1):min(((page)*npar),n)

	onoff{i} = on;

end

set_param(gcb,'MaskVisibilities',onoff);



