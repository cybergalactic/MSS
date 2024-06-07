function displayVehicleData(vehicleName, vehicleData, imageFile, figNo)
% Compatible with MATLAB and the free software GNU Octave (www.octave.org).
% displayVehicleData(vehicleName, vehicleData, imageFile, figNo) displays
% the vehicle main characteristics and an image of the vehicle.
%
% Inputs:
%   vehicleName - A string representing the name of the vehicle 
%                 (e.g., 'Remus 100 AUV')
%   vehicleData - A cell array of key-value pairs representing the vehicle's 
%                 characteristics (e.g., {'Length', '1.6 m', 
%                 'Diameter', '19 cm', 'Mass', '31.9 kg', ...})
%   imageFile   - A string representing the file name of the vehicle image 
%                 located in .../MSS/SIMULINK/figs/ (e.g., 'remus100.jpg')
%   figNo       - A number representing the figure number to use for the display
%
% Outputs:
%   None
%
% Author:     Thor I. Fossen
% Date:       2024-06-07

% Create the heading text
Heading = sprintf('%-30s\n%-30s\n', 'MSS Toolbox', vehicleName);

% Create the data text with proper alignment
numEntries = length(vehicleData) / 2;
formatSpec = repmat('%-30s : %s\n', 1, numEntries);
MSS_text = sprintf(formatSpec, vehicleData{:});

% Create a figure window
figure(figNo); figure(gcf);

% Create an axes for the heading
axes('Position', [0.1 0.8 0.8 0.1]);
text(0, 1, Heading, ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Courier', ...
    'HorizontalAlignment', 'left');
axis off;

% Create an axes for the text
axes('Position', [0.1 0.53 0.8 0.3]);
text(0, 1, MSS_text, ...
    'FontSize', 16, ...
    'FontName', 'Courier', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
axis off;

% Read and display the image
filePath = which(imageFile);
axes('Position', [0.1 0.05 0.8 0.4]);
imshow(imread(filePath));

end