function [pyPath, vers] = getPython()
    pyPath = "";
    vers = "";
    
    % Detect the operating system
    if ispc
        % Windows: try to use 'python' from the system PATH
        [status, cmdout] = system('where python');
        if status ~= 0 || isempty(cmdout)
            error('Python not found in the system PATH on Windows.');
        end
        pyList = strtrim(strsplit(cmdout, newline));
        pyPath = pyList{1};  % take the first result

        % Get Python version
        [~, cmdout] = system(['"' pyPath '" --version']);
        vers = strtrim(cmdout);
    else
        % Mac/Linux: try to use 'python3' from the system PATH
        [status, cmdout] = system('which python3');
        if status ~= 0 || isempty(cmdout)
            fprintf('Python3 not found in the system PATH on Mac/Linux.');
        end

        pyPath = strtrim(cmdout);
        [~, cmdout] = system('python3 --version');
        vers = strtrim(cmdout);
    end

    % Display the Python path being used
    %log = ['Setting Python to: ', pyPath];

    % Set pyenv
    % try
    %     %pyenv('Version', pyPath);
    %     % Show current Python environment
    % catch ME
    %     error('Failed to set pyenv: %s', ME.message);
    % end

end
