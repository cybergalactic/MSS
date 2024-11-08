function displayVesselStructure(vessel)
% displayVesselStructure prints all parameters in the vessel structure
fprintf('Vessel Parameters:\n');
fprintf('%-20s %-10s\n', 'Field', 'Value'); % Column headers
fprintf('%-20s %-10s\n', '--------------------', '----------'); % Separator line

% Loop through each field in the structure
fields = fieldnames(vessel);
for i = 1:numel(fields)
    fieldName = fields{i};
    fieldValue = vessel.(fieldName);
    
    % Check the type of the fieldValue for appropriate display
    if isnumeric(fieldValue) && isscalar(fieldValue)
        % Display numeric scalars
        fprintf('%-20s %-10g\n', fieldName, fieldValue);
    elseif ischar(fieldValue)
        % Display strings
        fprintf('%-20s %-10s\n', fieldName, fieldValue);
    elseif isnumeric(fieldValue) && isvector(fieldValue)
        % Display vectors as a single row
        fprintf('%-20s %-10s\n', fieldName, mat2str(fieldValue));
    elseif isnumeric(fieldValue) && ismatrix(fieldValue)
        % Display matrices with their full content
        fprintf('%-20s\n', fieldName);
        disp(fieldValue); % Print matrix content
    else
        % Placeholder for other complex types (cell arrays, structs)
        fprintf('%-20s %-10s\n', fieldName, '[Complex Value]');
    end
end

end