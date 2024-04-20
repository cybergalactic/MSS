function ControlFlag = controlMethod(methods)
% ControlFlag = controlMethod(methods) creates a GUI for selecting a control
% method from a list in MATLAB or a command-line menu in GNU Octave.
%
% Inputs:
%   methods - Cell array of strings, each representing a control method.
%
% Outputs:
%   ControlFlag - The index of the selected method.
%
% Compatibility:
%   Compatible with MATLAB and GNU Octave.
%
% Author:     Thor I. Fossen
% Date:       2021-04-25

if isoctave()  % Octave command line interface
    fprintf('\nSelect a control method from the list below:\n');
    for i = 1:length(methods)
        fprintf('%d: %s\n', i, methods{i});
    end
    validInput = false;
    while ~validInput
        idx = input('Enter your choice (number): ', 's');  % Read input as string
        idx = str2double(idx);  % Convert input string to double
        if ~isnan(idx) && idx == fix(idx) && idx >= 1 && idx <= length(methods)
            validInput = true;
            ControlFlag = idx;
        else
            fprintf('Invalid selection. Please enter a number between 1 and %d.\n', length(methods));
        end
    end
else % MATLAB GUI
    numMethods = length(methods);
    buttonHeight = 40;
    buttonSpacing = 10;
    totalHeight = numMethods * (buttonHeight + buttonSpacing) + buttonSpacing;
    figureWidth = 500;
    figureHeight = max(250, totalHeight);  % Ensure minimum height is 250

    % Create the figure window
    fig = figure('Name', 'Choose Control Method', 'Position', ...
        [200, 300, figureWidth, figureHeight], 'MenuBar', 'none', ...
        'NumberTitle', 'off', 'Resize', 'off', 'CloseRequestFcn', @closeDialog);

    fig.UserData = NaN;  % Initialize UserData to store ControlFlag

    % Create push buttons for each method
    for i = 1:numMethods
        posY = figureHeight - (i * (buttonHeight + buttonSpacing));
        uicontrol('Style', 'pushbutton', 'String', methods{i}, ...
            'FontSize', 13, 'Position', [10, posY, 450, 40], ...
            'Callback', {@buttonCallback, i});
    end

    uiwait(fig);                    % Wait for figure to close
    ControlFlag = fig.UserData;     % Retrieve ControlFlag
    delete(fig);                    % Delete figure to clean up
end

% Callback function for push buttons
function buttonCallback(~, ~, index)
    fig.UserData = index;       % Store selection in figure's UserData
    uiresume(fig);              % Resume execution of uiwait
end

% Close Request Function to handle user closing the figure window
function closeDialog(src, ~)
    uiresume(src); % Allow uiwait to return even if the window is closed
end

end

