function displayMfileHeader(filename)
    disp('--------------------------------------------------------------------');
    fprintf('MSS toolbox: %s:\n', filename)
    disp('--------------------------------------------------------------------');
    
    % Open the file
    fid = fopen(filename, 'r');
    
    % Check if the file opened successfully
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read and display lines starting with '%' until a non-comment line is encountered
    while ~feof(fid)
        line = fgetl(fid);
        
        % Disregard the 'function' line
        if startsWith(line, 'function')
            continue;
        end
        
        if startsWith(line, '%')
            % Check if the line is just a '%' or '% ' (space after %)
            if isscalar(line) || all(line(2:end) == ' ')
                disp(' ');  % Print an empty line
            else
                % Strip off the leading '%' and any following space
                headerLine = strtrim(line(2:end));
                disp(headerLine);
            end
        else
            break;
        end
    end
    
    % Close the file
    fclose(fid);

    disp('--------------------------------------------------------------------');
end