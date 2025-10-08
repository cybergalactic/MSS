function displayVehicleData(vehicleName, vehicleData, imageFile, figNo)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org).
% displayVehicleData(vehicleName, vehicleData, imageFile, figNo) displays
% the vehicle main characteristics and an image of the vehicle.
%
% Inputs:
%   vehicleName - A string representing the name of the vehicle 
%   vehicleData - A cell array of key-value pairs representing the vehicle's 
%                 characteristics
%   imageFile   - File name of the vehicle image (in path or MSS/SIMULINK/figs/)
%   figNo       - Figure number
%
% Author:     Thor I. Fossen
% Date:       2024-06-07

% Heading and table text
Heading = sprintf('%s\n%s\n', 'MSS Toolbox', vehicleName);
numEntries = numel(vehicleData)/2;
formatSpec = repmat('%-25s : %s\n', 1, numEntries);
MSS_text = sprintf(formatSpec, vehicleData{:});

% Figure window
figure(figNo); clf;
set(gcf, 'Name', vehicleName, 'NumberTitle', 'off', 'Color', 'w');
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [100 100 480 420]);   
set(gcf, 'Resize', 'off');

% Heading (top) 
axes('Position', [0.06 0.84 0.88 0.12]);
text(0, 0.95, Heading, ...
    'FontSize', 13, 'FontWeight', 'bold', ...
    'FontName', 'Courier', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
axis off;

% Vehicle info text (middle)
axes('Position', [0.06 0.56 0.88 0.25]);
text(0, 1, MSS_text, ...
    'FontSize', 11, 'FontName', 'Courier', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'Interpreter', 'none');
axis off;

% Image 
filePath = which(imageFile);
axes('Position', [0.15 0.20 0.7 0.33]);  
if isempty(filePath)
    text(0, 0.5, sprintf('Image not found:\n%s', imageFile), ...
        'FontSize', 11, 'Color', 'r', 'FontName', 'Courier', ...
        'Interpreter', 'none');
    axis off;
else
    img = imread(filePath);
    imshow(img, 'InitialMagnification', 'fit');
end

end