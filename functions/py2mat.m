function out = py2mat(obj)
%PY2MAT Converts Python objects returned to MATLAB into native MATLAB types.
%
%   out = PY2MAT(obj)
%
%   Supported conversions:
%       py.dict        -> struct
%       py.list/tuple  -> cell array
%       py.str         -> char
%       py.int/float   -> double
%       py.bool        -> logical
%       py.None        -> []
%       numpy.ndarray  -> double array (if numeric)
%       Otherwise leaves object as-is

    % ----- Handle Python None -----
    if strcmp(class(obj), 'py.NoneType')
        out = [];
        return;
    end

    % ----- Handle Python dict -----
    if isa(obj, 'py.dict')
        % Convert dict -> struct
        pykeys = py.list(obj.keys());      % convert dict_keys → list Python
        keys = cell(pykeys);               % convert list → MATLAB cell array

        out = struct();
        for i = 1:numel(keys)
            key = string(keys{i});
            value = obj{keys{i}};          % now safe indexing
            out.(key) = py2mat(value);     % recursive conversion
        end
        return;
    end

    % ----- Handle Python list or tuple -----
    if isa(obj, 'py.list') || isa(obj, 'py.tuple')
        out = cellfun(@py2mat, cell(obj), 'UniformOutput', false);
        return;
    end

    % ----- Handle strings -----
    if isa(obj, 'py.str')
        out = char(obj);
        return;
    end

    % ----- Handle Python numeric types -----
    if isa(obj, 'py.int') || isa(obj, 'py.float')
        out = double(obj);
        return;
    end

    % ----- Handle Python bool -----
    if isa(obj, 'py.bool')
        out = logical(obj);
        return;
    end

    % ----- Handle numpy ndarray -----
    if contains(class(obj), 'numpy.ndarray')
        try
            out = double(obj.tolist());   % convert to MATLAB double
        catch
            out = cellfun(@py2mat, cell(obj.tolist()), 'UniformOutput', false);
        end
        return;
    end

    % ----- Handle other Python objects (fallback) -----
    try
        % try to convert to double
        out = double(obj);
    catch
        try
            % try to convert to char
            out = char(obj);
        catch
            % last resort: keep as Python object
            out = obj;
        end
    end
end
