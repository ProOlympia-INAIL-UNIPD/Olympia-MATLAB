function S = renameFieldsRecursive(S, oldStr, newStr)
%RENAMEFIELDSRECURSIVE  Replace substring in all field names of a structure (recursively)
%
%   S = renameFieldsRecursive(S, oldStr, newStr)
%
%   Inputs:
%     S       - struttura arbitrariamente annidata
%     oldStr  - sottostringa da sostituire nei nomi dei campi
%     newStr  - nuova sottostringa da inserire
%
%   Output:
%     S       - struttura con i campi rinominati ricorsivamente

    if ~isstruct(S)
        return;
    end

    oldFields = fieldnames(S);
    for i = 1:numel(oldFields)
        field = oldFields{i};
        newField = strrep(field, oldStr, newStr);

        % Rinomina se necessario
        if ~strcmp(field, newField)
            [S.(newField)] = S.(field);
            S = rmfield(S, field);
        end

        % Ricorsione per eventuali sotto-strutture
        if isstruct(S.(newField))
            if numel(S.(newField)) > 1
                for k = 1:numel(S.(newField))
                    S.(newField)(k) = renameFieldsRecursive(S.(newField)(k), oldStr, newStr);
                end
            else
                S.(newField) = renameFieldsRecursive(S.(newField), oldStr, newStr);
            end
        end
    end
end