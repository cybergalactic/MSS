% GDF_EDIT  script for generation of a new Wamit GDF file fixing the problem
% that the waterline is not a closed curve (edit the file manually)
%
% Authors: Thor I. Fossen and �yvind Smogeli
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

filename = 'myship'  % *.gdf

eps_z    = 0.05;  % resolution for detecting points on the waterline z=0
eps_dist = 0.1;   % resolution for detecting neighboring points

%--------------------------------------------------------------------------
% plot original GDF data in Figures 1 and 2
%--------------------------------------------------------------------------
gdf_data = plot_wamitgdf(filename,'b',1,1);

%--------------------------------------------------------------------------
% Find points on the surface, merge close points and store these in
% the tabke gdf_new
%--------------------------------------------------------------------------
S = [];
Npanels = max(size(gdf_data));

for n = 1:Npanels
    S = [S ; gdf_data(n,1:3); gdf_data(n,4:6); gdf_data(n,7:9); gdf_data(n,10:12)];
end

Npoints = max(size(S));
for m = 1:Npoints
    if abs(S(m,3)) < eps_z         
        for k = 1:Npoints      
            dist = norm(S(k,:)-S(m,:));         
            if dist < eps_dist        
                S(k,:) = S(m,:);               
            end              
        end
    end    
end

k = 0;
for n = 1:Npanels
    temp = [];
    for m = 1:4
        k = k + 1;
        temp = [temp S(k,:)];
    end
    gdf_new(n,:) = temp;
end

%-------------------------------------------------------------------------
% save modified data as 'filname'new.gdf
%-------------------------------------------------------------------------
fid1 = fopen([filename '.gdf']);            % open Wamit file
fid2 = fopen([filename 'new.gdf'],'w');

txt1 = fgetl(fid1);
txt2 = fgetl(fid1);
txt3 = fgetl(fid1);
txt4 = fgetl(fid1);

fprintf(fid2,'%s\n',' MC generated Wamit GDF file using ShipX data');
fprintf(fid2,'%s\n',txt2);
fprintf(fid2,'%s\n',txt3);
fprintf(fid2,'%s\n',txt4);

for i=1:Npanels
    txt = num2str(gdf_new(i,:));
    fprintf(fid2,' %s\n',txt);
end

fclose(fid1);
fclose(fid2);

%-------------------------------------------------------------------------
% plot new data in Figures 3 and 4
%-------------------------------------------------------------------------
plot_wamitgdf([filename 'new'],'g',3,1);


