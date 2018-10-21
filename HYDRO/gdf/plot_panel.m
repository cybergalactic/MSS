function plot_panel(no,gdf_data,colorcode)
% PLOT_PANEL (MSS HYDRO)
%
% Plot Wamit panels from *.gdf file, see plot_wamitgdf.m
%
% >> plot_panel(n,gdf_data,colorcode)
%
% -------------------------------------------------------------------------
% Inputs:
%    no:           line number to be plotted, that is gdf_data(no,:)
%    colorcode:    'r','g','b','c','m','y','w','k' or colormaps [0 0 1] etc.
%    gdf_data(:,:) table dimension N x 12 of Wamit panel data
%
% -------------------------------------------------------------------------
% Author:    Thor I. Fossen and Ø. N. smogeli
% Date:      2005-09-01
% Revisions: 
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

panel = gdf_data(no,:);
x{no} = [panel(1), panel(4), panel(7), panel(10)];
y{no} = [panel(2), panel(5), panel(8), panel(11)];
z{no} = [panel(3), panel(6), panel(9), panel(12)];

patch(x{no},y{no},z{no},colorcode);
