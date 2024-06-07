function ControlFlag = controlMethod(methods)
% ControlFlag is compatible with MATLAB and GNU Octave (www.octave.org). 
% ControlFlag = controlMethod(methods) creates a GUI for selecting a control
% method from a list.
%
% Inputs:
%   methods - Cell array of strings, each representing a control method.
%
% Outputs:
%   ControlFlag - The index of the selected method.
%
% Compatibility:
%   Compatible with MATLAB and GNU Octave 9.1.0 and later.
%
% Author:     Thor I. Fossen
% Date:       2024-04-25
%   2024-05-09 : Same GUI for Matlab and GNU Octave. 

numMethods = length(methods);
buttonHeight = 40;
buttonSpacing = 10;
totalHeight = numMethods * (buttonHeight + buttonSpacing) + buttonSpacing;
figureWidth = 500;
figureHeight = max(150, totalHeight);  % Ensure minimum height is 150

% Create the figure window
f = figure('Name', 'Choose Control Method', 'Position', ...
    [200, 300, figureWidth, figureHeight], 'MenuBar', 'none', ...
    'NumberTitle', 'off', 'Resize', 'off', 'CloseRequestFcn', ...
    @(src, event) closeDialog(src, event));

set(f, 'UserData', NaN);  % Initialize UserData to store ControlFlag

% Nested Callback Function
    function buttonCallback(src, ~)
        % Check if f is a figure handle
        if ishandle(f) && strcmp(get(f, 'Type'), 'figure')  
            set(f, 'UserData', get(src, 'UserData'));      % Get the index
            uiresume(f);
        else
            disp('Error: f is not a valid figure handle.');
        end
    end

% Create push buttons for each method
for i = 1:numMethods
    posY = figureHeight - (i * (buttonHeight + buttonSpacing));
    uicontrol('Style', 'pushbutton', 'String', methods{i}, ...
        'UserData', i, 'FontSize', 13, 'Position', [10, posY, 450, 40], ...
        'Callback', @buttonCallback);
end

drawnow;  % Forces a redraw of the GUI to update its appearance
uiwait(f);                          % Wait for figure to close
ControlFlag = get(f, 'UserData');   % Retrieve ControlFlag using get
delete(f);                          % Delete figure to clean up

end

% Close Request Function to handle user closing the figure window
function closeDialog(src, ~)
    uiresume(src); % Allow uiwait to return even if the window is closed
end

