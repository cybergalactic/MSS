function gdf2body(vessel_name,code)
% GDF2BODY (MSS HYDRO)
%
% gdf2body(vessel_name,code) transforms the GDF file from h-frame origin to the
% b-frame origin. It is also possible to modify the acceleration of gravity in
% the GDF file. All data in WAMIT are computed in the b-frame.
%
% Example:
%
% >> gdf2body('supply','seaway') updates the GDF file supply.gdf and stores 
%                                theresults in supplyCG.gdf
%
% Function input: vessel_name    string describing the GDF-file
%                 code           string: 'veres' or 'seaway'
%
% -------------------------------------------------------------------------
% 
% User inputs:    g         Acceleration of gravity (m/s^2)
%                 xg,yg,zg  Body-fixed CG with respect to the global coordinate
%                           system (Lpp/2,0,WL). Axes: forward, port, upwards
%
% Output:         vesselCG  updated GDF file with new CG and g
% -------------------------------------------------------------------------
% Author:    Thor I. Fossen
% Date:      2005-09-17
%
% Revisions: 2005-10-20 (TIF): extended to read octopus seaway files
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

GDFfile = [vessel_name '.gdf'];

% Update acceleration of gravity 
g = input('Acceleration of gravity - default 9.81 m/s^2: ');
if ~exist(num2str(g))
    g = 9.81;
end

% Update CG 
disp(' ')
disp('CG w.r.t. to the global coordinate system: [Lpp/2 0 WL]');
CG = input('with axes forward, port, and upwards - default 0 0 0: ','s');

if length(CG)==0
    CG = [0 0 0];
else
    CG = str2num(CG);
end

% open GDF-file
if ~exist(GDFfile)
    error('Error: WAMIT GDF-file does not exist')
    return
end

fid1 = fopen(GDFfile);

% read GDF file
if strcmp(code,'veres')

    % 4 lines before data
    txt1 = fgetl(fid1);
    txt2 = fgetl(fid1);
    data = str2num(txt2);  % update acceleration of gravity
    txt2 = [' ' num2str([data(1) g])];
    txt3 = fgetl(fid1);
    txt4 = fgetl(fid1);
    Npanels = str2num(txt4);

    for i=1:Npanels
        panel_data(i,:) = str2num(fgetl(fid1));
        % update GDF file with new CG
        panel_data(i,:) = panel_data(i,:)  - [CG CG CG CG];
    end

    fclose(fid1);

    % save modified file

    GDFfile2 = [vessel_name 'CG.gdf'];
    fid2 = fopen(GDFfile2,'w');
    fprintf(fid2,'%s\n',txt1');
    fprintf(fid2,'%s\n',txt2);
    fprintf(fid2,'%s\n',txt3);
    fprintf(fid2,'%s\n',txt4);

    for i=1:Npanels
        fprintf(fid2,' %f %f %f %f %f %f %f %f %f %f %f %f \n',panel_data(i,:));
    end

elseif strcmp(code,'seaway')   

    % 7 lines before data
    txt1 = fgetl(fid1);
    txt2 = fgetl(fid1);
    txt3 = fgetl(fid1);
    txt4 = fgetl(fid1);    
    txt5 = fgetl(fid1);      
    data = str2num(txt5);  % update acceleration of gravity
    txt5 = [' ' num2str([data(1) g])];
    txt6 = fgetl(fid1);
    txt7 = fgetl(fid1);
    Npanels = str2num(txt7);

    for i=1:4*Npanels  % 4 lines = 1 panel
        panel_data(i,:) = str2num(fgetl(fid1));
    end

    % compute x-distance from saway origin to Lpp/2
    Loa = max(panel_data(:,1))-min(panel_data(:,1));
    disp(sprintf('Length overall computed from GDF data: L_oa = %2.2f',Loa))
    Lpp  = input('Length between perpendiculars        : L_pp = ');
    x_globalCO = Lpp/2 - min(panel_data(:,1));
    
    for i=1:4*Npanels  % update GDF file with new CG
        panel_data(i,:) = panel_data(i,:) -[x_globalCO 0 0] - CG;
    end

    fclose(fid1);

    % save modified file

    GDFfile2 = [vessel_name 'CG.gdf'];
    fid2 = fopen(GDFfile2,'w');
    fprintf(fid2,'%s\n',txt1');
    fprintf(fid2,'%s\n',txt2);
    fprintf(fid2,'%s\n',txt3);
    fprintf(fid2,'%s\n',txt4);
    fprintf(fid2,'%s\n',txt5);
    fprintf(fid2,'%s\n',txt6);
    fprintf(fid2,'%s\n',txt7);    

    for i=1:4*Npanels
        fprintf(fid2,' %f %f %f  \n',panel_data(i,:));
    end
    
else
    
    disp('Error: wrong input argument')
    return
    
end

% plot new data
figure(1),clf(1)
plot_wamitgdf([vessel_name 'CG'],'blue',1);

zz = CG(3);
xx = 60;
yy = 20;

panel  = [-xx -yy zz xx -yy zz xx yy zz -xx yy zz];
x = [panel(1), panel(4), panel(7), panel(10)];
y = [panel(2), panel(5), panel(8), panel(11)];
z = [panel(3), panel(6), panel(9), panel(12)];

patch(x,y,z,'green');

fclose(fid2);

% display results
disp(' ')
disp('The GDF file has been updated.');
disp(['New CG is: [' num2str(CG) '] w.r.t. [Lpp/2 0 WL]'])
disp(['New acceleration of gravity is: ' num2str(g) '(m/s^2)'])
disp(' ')






