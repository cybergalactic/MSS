function plot_wamitidf(filename,colorcode,figno)
% PLOT_WAMITIDF (MSS HYDRO)
%
% Plot Wamit geometry files (*.gdf) of ships and rigs in 3D
%
% Examples:
%
% >> plot_wamitidf(filename,colorcode,figno)
%
% >> plot_wamitidf(filename);  plot filname.gdf in Figure 1 with default color
%
% -------------------------------------------------------------------------
% Inputs:
%    filename:    Wamit IDF file without extension
%    colorcode:   'r','g','b','c','m','y','w','k' or colormaps [0 0 1] etc.
%    figno:       figure number
%
% -------------------------------------------------------------------------
% Author:    Thor I. Fossen and Ø. N. smogeli
% Date:      2005-09-06
% Revisions: 
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

if ~exist('colorcode')
    colorcode = 'b';
end

if ~exist('figno')
    figno = 1;
end

% -------------------------------------------------------------------------
% Read Wamit panels
% -------------------------------------------------------------------------
fid1 = fopen([filename '.idf']);   % open Wamit file

txt1 = fgetl(fid1);
txt2 = fgetl(fid1);
txt3 = fgetl(fid1);
txt4 = fgetl(fid1);

Npanels = str2num(txt4);

i = 1;
    
while feof(fid1) == 0,
	
	for j = 1:4

		txt       = char(fgetl(fid1));
		numbers   = str2num(txt);
		idf_data(i,((j-1)*3+1):j*3) = numbers;

	end

	i = i+1;

end

fclose(fid1);                    % close Wamit file

% -------------------------------------------------------------------------
% Plot Wamit panels in 3D
% -------------------------------------------------------------------------

figure(figno)

% Loop over panels
for n = 1:Npanels
	
	% Vector with panel coordinates
	panel = idf_data(n,:);
	x{n} = [panel(1), panel(4), panel(7), panel(10)];
    y{n} = [panel(2), panel(5), panel(8), panel(11)];
    z{n} = [panel(3), panel(6), panel(9), panel(12)];

    % plot vessel
    if ~exist('half') | half ~= 1
        patch(x{n},y{n},z{n},colorcode);
        patch(x{n},-y{n},z{n},colorcode);
    else
        patch(x{n},y{n},z{n},colorcode);
    end
end

axis('equal')
grid on
xlabel('X-axis (m)')
ylabel('Y-axis (m)')
zlabel('Z-axis (m)')
title(strcat('3D Visualization of the auto-generated Wamit file:',[' ',filename '.idf']))
view(45,16)
