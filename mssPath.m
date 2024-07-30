% MSS Path Update Script
% This script checks for the presence of the 'MSS' directory on MATLAB's 
% path. If found, it removes all existing instances of this path and its 
% subdirectories and then reinstalls them to ensure the path includes all 
% current subdirectories. This is useful for situations where the directory 
% structure of 'MSS' may change due to updates or modifications in the 
% file system.
%
% Usage:
% - Run this script to refresh the MATLAB path entries related to the
%   'MSS' directory.
% - It handles multiple instances of 'MSS' directories by using the
%   first found.
% - It temporarily suppresses and then restores warnings related to 
%   non-existent path entries to prevent clutter in the command window.
%
% Author:    Thor I. Fossen
% Date:      2024-04-26
% Revisions:

if exist('OCTAVE_VERSION', 'builtin')
    disp(['Add the MSS path with subfolders by starting' ...
        ' "octave --gui", then click, edit - set path']);
    return;
end

% MATLAB specific warning suppression
warning('off', 'MATLAB:rmpath:DirNotFound');

% Locate directories named 'MSS' using the MATLAB what function
mssInfo = what('MSS');

% Evaluate the presence and number of 'MSS' directories found
if isempty(mssInfo)
    error(['MSS directory not found, please add MSS path with' ...
        ' subfolders from the menu.']);
elseif length(mssInfo) > 1
    % Handle multiple MSS directories found by notifying the user
    % and using the first one
    warning(['%d MSS directories found. Using the first one' ...
        ' found at: %s'], length(mssInfo), char(mssInfo(1).path));
    % Use the first directory found, ensure it's char
    basePath = char(mssInfo(1).path);  
else
    % Confirm when one MSS path is found
    basePath = char(mssInfo.path);  % Ensure the path is char
    fprintf('MSS directory found at: %s\n', basePath);
end

% Remove old MSS paths to clean up any outdated links
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath(basePath));
warning('on', 'MATLAB:rmpath:DirNotFound');

% Re-add the MSS directory and all its subdirectories to MATLAB's path
addpath(genpath(basePath));

% Save the updated path
savepath;
fprintf('MATLAB path updated and saved successfully in MATLAB.\n');