function pyPath = setDefaultPython()
    % Detect the operating system
    if ispc
        % Windows: try to use 'python' from the system PATH
        [status, cmdout] = system('where python');
        if status ~= 0
            error('Python not found in the system PATH on Windows.');
        end
        pyList = strtrim(strsplit(cmdout, newline));  % take the first result
        pyPath = pyList{1};
    else
        % Mac/Linux: try to use 'python3' from the system PATH
        [status, cmdout] = system('which python3');
        if status ~= 0
            error('Python3 not found in the system PATH on Mac/Linux.');
        end
        pyPath = strtrim(cmdout);
    end

    % Display the Python path being used
    fprintf('Setting Python to: %s\n', pyPath);

    % Set pyenv
    try
        pyenv('Version', pyPath);
    catch ME
        error('Failed to set pyenv: %s', ME.message);
    end

    % Show current Python environment
    %disp(pyenv);
end
