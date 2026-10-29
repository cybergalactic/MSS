function gdf_data = plot_wamitgdf(filename,colorcode,figno,waterline,half)
% PLOT_WAMITGDF (MSS Hydro)
%
% Plot Wamit low-order geometry files (*.gdf) of ships and rigs in 3D.
% GDF-file format:
%
% header
% ULEN GRAV
% ISX ISY
% NPAN
% X1(1) Y1(1) Z1(1) X2(1) Y2(1) Z2(1) X3(1) Y3(1) Z3(1) X4(1) Y4(1) Z4(1)
% X1(2) Y1(2) Z1(2) X2(2) Y2(2) Z2(2) X3(2) Y3(2) Z3(2) X4(2) Y4(2) Z4(2)
% .
% . . . . . . . . . . . . . . . . . . . . . . . X4(NPAN) Y4(NPAN) Z4(NPAN)
%
% Examples:
%
% >> plot_wamitgdf(filename,colorcode,figno,waterline,half);
%
% >> plot_wamitgdf(filename);  plot filname.gdf in Figure 1 with default color
% >> gdf_data = plot_wamitgdf(filename);  returns data in table gdf_data
%
% Use the following command to check/debug the free surface points needed to
% remove irregular frequencies in Wamit (IRR = 3, see the Wamit manual)
%
% >> plot_wamitgdf(filename,colorcode,figno,1);
% -------------------------------------------------------------------------
% Inputs:
%    filename:    Low-order Wamit GDF file without extension
%    colorcode:   'r','g','b','c','m','y','w','k' or colormaps [0 0 1] etc.
%    figno:       figure number
%    waterline:   optional flag used to plot free surface points
%    half:        optional, use 1 to plot GDF data (no mirror)
%                           for Veres data. Octopus data always has mirror
%
% Output:
%    gdf_data:    Table dimension N x 12 of Wamit panels where N is the 
%                 number of panels given by 3x4=12 points on each line
%
% -------------------------------------------------------------------------
% Author: Thor I. Fossen and Ø. N. smogeli
% Date:   2005-09-01 
%
% Revisions: 2005-10-20 (TIF): handles both WAMIT GDF-data formats, e.g.:
%                              1 point on each line
%                              4 points on each line
%            2008-01-12 (TIF): imrpoved handling of text strings
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

clear gdf_data x y z

if ~exist('colorcode')
    colorcode = [0.3 0.5 0.9];
end

if ~exist('figno')
    figno = 1;
end

if ~exist('waterline')
    waterline = 0;
end

% -------------------------------------------------------------------------
% Read Wamit panels
% -------------------------------------------------------------------------
fid1 = fopen([filename '.gdf']);   % open Wamit file

header = fgetl(fid1);
disp(header)
txt = fgetl(fid1);
if isempty(findstr('ULEN',txt)) == 1
    gdata = str2num(txt);
else
    Nt = findstr('ULEN',txt);
    gdata = str2num(txt(1:Nt-1));
end
disp(['ULEN = ' num2str(gdata(1))]);
disp(['GRAV = ' num2str(gdata(2))]);

txt = fgetl(fid1);
if isempty(findstr('ISX',txt)) == 1
    symmetry = str2num(txt);
else
    Nt = findstr('ISX',txt);
    symmetry = str2num(txt(1:Nt-1));
end
disp(['ISX = ' num2str(symmetry(1))]);
disp(['ISY = ' num2str(symmetry(2))]);
   
txt = fgetl(fid1);
if isempty(findstr('NPAN',txt)) == 1
    Npanels = str2num(txt);
else
    Nt = findstr('NPAN',txt);
    Npanels = str2num(txt(1:Nt-1));
end
disp(['NPAN = ' num2str(Npanels)]);

% check if first data line has 3 or 12 elements
txt      = fgetl(fid1);
numbers1 = str2num(txt);
disp(['GDF file format: ' num2str(length(numbers1)) ' numbers on each line']);
disp('Rading data, please wait....')

% store first line in gdf_data
if length(numbers1) == 12,
    format = 1;    
    gdf_data(1,1:12) = numbers1;       
else
    length(numbers1) == 3;
    format = 2;
    txt      = fgetl(fid1);
    numbers2 = str2num(txt);
    txt      = fgetl(fid1);
    numbers3 = str2num(txt);
    txt      = fgetl(fid1);
    numbers4 = str2num(txt);
    gdf_data(1,1:12) = [numbers1 numbers2 numbers3 numbers4];     
end

% lines 2:N
i = 2;

while feof(fid1) == 0,

    txt        = char(fgetl(fid1));
    numbers1   = str2num(txt); 
    
    if format == 1
        gdf_data(i,1:12) = numbers1;
    else  % format == 2
        txt      = fgetl(fid1);
        numbers2 = str2num(txt);
        txt      = fgetl(fid1);
        numbers3 = str2num(txt);
        txt      = fgetl(fid1);
        numbers4 = str2num(txt);
        gdf_data(i,1:12) = [numbers1 numbers2 numbers3 numbers4];
    end
    
    i = i+1;
    
end

fclose(fid1);                    % close Wamit file

% -------------------------------------------------------------------------
% Plot Wamit panels in 3D
% -------------------------------------------------------------------------
figure(figno)
clf(figno)
figure(gcf)

% Loop over panels
for n = 1:Npanels
	
	% Vector with panel coordinates
	panel = gdf_data(n,:);
	x{n} = [panel(1), panel(4), panel(7), panel(10)];
    y{n} = [panel(2), panel(5), panel(8), panel(11)];
    z{n} = [panel(3), panel(6), panel(9), panel(12)];

    % plot vessel
    if ~exist('half') | half ~= 1
        if symmetry(2) == 1, % Seaway generates half the ship                  
            patch(x{n},y{n},z{n},colorcode);
            patch(x{n},-y{n},z{n},colorcode);
        else % Octopus Modeller generates both ship sides
            patch(x{n},y{n},z{n},colorcode);
        end
    else % if half = 1, half the ship is plotted (for Veres only)
        patch(x{n},y{n},z{n},colorcode);
    end
end

axis('equal')
grid
xlabel('X-axis (m)')
ylabel('Y-axis (m)')
zlabel('Z-axis (m)')
title(strcat('3D Visualization of the Wamit file:',[' ',filename '.gdf']))
view(45,16)

% -------------------------------------------------------------------------
% Plot points on waterline - used to check free surface
% -------------------------------------------------------------------------
if waterline == 1,
    
    z_eps = 0.001;  % epsilon for finding waterline

    x_w = []; y_w = []; z_w = [];

    figure(figno+1)
    clf(figno+1)
    figure(gcf)
    
    for n = 1:Npanels

        index = find(abs(z{n}) < z_eps);

        if ~isempty(index)

            x_w = [x_w, x{n}(index)];
            y_w = [y_w, y{n}(index)];
            z_w = [z_w, z{n}(index)];

        end

    end

    [x_w,index] = sort(x_w);
    y_w = y_w(index);
    z_w = z_w(index);
    i = 0;

    disp(' ')
    disp('Waterline defined by only one point at:')
    clear x_err y_err z_err

    for n = 1:length(x_w);

        if n ~= 1 & n ~= length(x_w)

            test = find(x_w(n) == x_w);

            if size(test) < 2

                disp(sprintf('  x = %0.4f, y = %0.4f, z = %0.4f',x_w(n), y_w(n), z_w(n)))
                i = i + 1;
                x_err(i) = x_w(n);
                y_err(i) = y_w(n);
                z_err(i) = z_w(n);

            end

        end

    end

    plot(x_w,y_w,'.b')
    hold on
    plot(x_w,-y_w,'.b')

    if i > 0
        plot(x_err,y_err,'or')
    end

    axis equal
    grid
    hold off

end
