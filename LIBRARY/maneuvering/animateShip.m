function animateShip(xPath, yPath, shipSize, lineColor, figNo)
% animateShip is compatibel with MATLAB but not supported in GNU Octave.
% The function animateShip(xPath, yPath, shipSize, lineColor, figNo) 
% animates a viking ship moving along a specified path on a North-East 
% (y-x) plot. The function takes a path defined by North (yPath) and East 
% (xPath) coordinates,animates a viking ship image moving along this path, 
% and scales the ship size based on the specified proportion of the plot 
% width. The path and ship are plotted on a figure identified by figNo, 
% allowing for visualization of navigational paths. It is recommended to 
% remove the persistent variables by adding 'clear animateShip' on top of 
% your script.
%
% Inputs:
%   xPath     - x coordinates of the path (1xN vector)
%   yPath     - y coordinates of the path (1xN vector)
%   shipSize  - Proportion of the plot width that the ship's width should 
%               occupy (scalar)
%   lineColor - Color of the path line (string or RGB vector)
%   figNo     - Figure number for plotting (scalar)
%
% Example usage:
%   animateShip([0 1 2 3], [0 1 4 9], 0.1, 'b', 1)
%   This will animate a ship along a specified path in figure 1, with the 
%   path colored blue and the ship's width being 5% of the plot's width.
%
% Dependencies:
%   Requires the image file 'viking.png' located in the MSS/GNC/utils 
%   directory. 
%
% Author:    Thor I. Fossen
% Date:      2024-03-27
% Revisions: 
%   2024-05-09 : Use which('viking.png') to find the path in Matlab and Octave.
%
% Note: Downsampling is applied if the path contains more than 200 points 
% to improve performance. The ship is automatically scaled and oriented to 
% follow the path direction.

persistent shipImgRotated shipAlphaRotated shipHandle

if isempty(shipImgRotated)

    filePath = which('viking.png');
    if isempty(filePath)
        error('The file viking.png not found on the MSS path.');
    else
        [shipImg, ~, shipAlpha] = imread(filePath);  % MSS/GNC/utils/..
        shipImgRotated = rot90(shipImg, 2); % Rotate if necessary
        if ~isempty(shipAlpha)
            shipAlphaRotated = rot90(shipAlpha, 2);
        else
            shipAlphaRotated = [];
        end
    end
end

% Determine downsampling if needed
N = length(xPath);
if N > 200
    Nsamples = round(N / 200);
    xPath = xPath(1:Nsamples:end);
    yPath = yPath(1:Nsamples:end);
end

figure(figNo);
clf; % Clear figure once at the start

% Plot the path once at the start
plot(yPath, xPath, lineColor, 'LineWidth', 2);
hold on;
grid on;

xlabel('East', 'FontSize', 14);
ylabel('North', 'FontSize', 14);
title('North-East positions', 'FontSize', 14);
axis equal;

% Initialize the ship image handle if it doesn't exist or is invalid
if isempty(shipHandle) || ~isgraphics(shipHandle, 'image')
    shipHandle = image('CData', shipImgRotated, 'AlphaData', shipAlphaRotated);
end

% Calculate scaling for the ship based on plot dimensions
plotWidth = max(yPath) - min(yPath);
shipWidth = plotWidth * shipSize;
aspectRatio = size(shipImgRotated, 2) / size(shipImgRotated, 1);
shipHeight = shipWidth / aspectRatio;

for i = 1:length(xPath)
    % Update ship position
    shipXCenter = yPath(i);
    shipYCenter = xPath(i);
    shipXData = [shipXCenter - shipWidth/2, shipXCenter + shipWidth/2];
    shipYData = [shipYCenter - shipHeight/2, shipYCenter + shipHeight/2];

    set(shipHandle, 'XData', shipXData, 'YData', shipYData);

    drawnow;
end

hold off;

end
