% GDF_WATERLINE script for generation of a new Wamit low-order GDF file
% where points above the waterline are removed.
%
% Author: Thor I. Fossen
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

filename = 'tanker'  % *.gdf
draft  = 7;          % desired  draft w.r.t. base line, remove panels above this

%--------------------------------------------------------------------------
% plot original GDF data in Figures 1
%--------------------------------------------------------------------------
%gdf_data = plot_wamitgdf(filename,'b',1);
[Npanels,dummy] = size(gdf_data);

% shift z-axis to new draft (GDF file should have z=0 in waterline)
z_base = min(gdf_data(:,3));   % baseline z-coordinate in GDF file

for i = 1:Npanels
    gdf_data(i,3)  = gdf_data(i,3)  - z_base - draft;
    gdf_data(i,6)  = gdf_data(i,6)  - z_base - draft;
    gdf_data(i,9)  = gdf_data(i,9)  - z_base - draft;    
    gdf_data(i,12) = gdf_data(i,12) - z_base - draft;   
end
%--------------------------------------------------------------------------
%header
%ULEN GRAV
%ISX ISY
%NPAN
%X1(1) Y1(1) Z1(1) X2(1) Y2(1) Z2(1) X3(1) Y3(1) Z3(1) X4(1) Y4(1) Z4(1)
%X1(2) Y1(2) Z1(2) X2(2) Y2(2) Z2(2) X3(2) Y3(2) Z3(2) X4(2) Y4(2) Z4(2)
%.
%. . . . . . . . . . . . . . . . . . . . . . . X4(NPAN) Y4(NPAN) Z4(NPAN)
%
% The vertices must be numbered in the counter-clockwise direction when the
% panel is viewed from the fluid domain.
%--------------------------------------------------------------------------

% check to see if numbering av vertices is OK - NOT IMPLEMENTED
disp('Warning: GDF data must be numbered counter-clockwise, this is not tested')

%--------------------------------------------------------------------------
% project panels above the waterline on the waterline
%--------------------------------------------------------------------------
gdf_new = [];
k = 1;

% waterline plane, z positive upwards
WL = createPlane([0 0 0],[0 0 1]);

for i = 1:Npanels

    p{1} = gdf_data(i,1:3);
    p{2} = gdf_data(i,4:6); 
    p{3} = gdf_data(i,7:9); 
    p{4} = gdf_data(i,10:12);
    
    % criterion for waterline test
    criterion = ([p{1}(3),p{2}(3),p{3}(3),p{4}(3)] < 0);
    
    % all vertices under the waterline
    if sum(criterion) == 4 
      
        gdf_new(k,:) = [p{1} p{2} p{3} p{4}];
        k = k+1;     

    % 1, 2 and 3 points under waterline    
    elseif sum(criterion) == 1 | sum(criterion) == 2 | sum(criterion) == 3
        if sum(cross(p{1},p{2})) == 0           % 3 vertices / triangles
            POLY1 = [p{1}; p{3}; p{4}];
        elseif sum(cross(p{1},p{3})) == 0
            POLY1 = [p{1}; p{2}; p{4}];
        elseif sum(cross(p{1},p{4})) == 0
            POLY1 = [p{1}; p{2}; p{3}];
        elseif sum(cross(p{2},p{3})) == 0
            POLY1 = [p{1}; p{2}; p{4}]; 
        elseif sum(cross(p{2},p{4})) == 0
            POLY1 = [p{1}; p{2}; p{3}];
        elseif sum(cross(p{3},p{4})) == 0
            POLY1 = [p{1}; p{2}; p{4}];            
        else 
            POLY1 = [p{1}; p{2}; p{3}; p{4}];    % 4 vertices
        end
        
        POLY2 = clipConvexPolygon3dPlane(POLY1, WL);
        [mx,my] = size(POLY1);

        % add 4th point for triangles, THIS IS NOT ROBUST
        if mx == 3  
            POLY2(4,:) = POLY2(3,:);
        end
        
        gdf_new(k,:) = [POLY2(1,:) POLY2(2,:) POLY2(3,:) POLY2(4,:)];
        k = k+1;    
        
    end
           
end

[Npanels_new, dummy] = size(gdf_new);

%-------------------------------------------------------------------------
% save modified data as 'filname'new.gdf
%-------------------------------------------------------------------------
fid1 = fopen([filename '.gdf']);            % open Wamit file
fid2 = fopen([filename 'new.gdf'],'w');

txt1 = fgetl(fid1);
txt2 = fgetl(fid1);
txt3 = fgetl(fid1);
txt4 = fgetl(fid1);

fprintf(fid2,'%s\n',['Modified GDF file with draft: ' num2str(draft) ' m' ]);
fprintf(fid2,'%s\n',txt2);                  % old data
fprintf(fid2,'%s\n',txt3);                  % old data
fprintf(fid2,'%s\n',num2str(Npanels_new));  % new panel number

for i=1:Npanels_new
    fprintf(fid2,' %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',gdf_new(i,:));
end

fclose(fid1);
fclose(fid2);

%-------------------------------------------------------------------------
% plot new data in Figures 2
%-------------------------------------------------------------------------
plot_wamitgdf([filename 'new'],'g',2);


