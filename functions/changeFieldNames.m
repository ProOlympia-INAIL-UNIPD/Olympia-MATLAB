function S = changeFieldNames(S, toremove, new)
if ~isstruct(S)
    return;
end

if nargin < 3
    new = "";
end

oldFields = fieldnames(S);
for i = 1:numel(oldFields)
    field = oldFields{i};
    newField = replace(field, toremove, new);
    
    if ~strcmp(field, newField)
        [S.(newField)] = S.(field);
        S = rmfield(S, field);
    end
end

end