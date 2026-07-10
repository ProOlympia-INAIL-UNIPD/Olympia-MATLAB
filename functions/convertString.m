function db_out = convertString(db_in)
    columnNames = db_in.Properties.VariableNames;
    
    db_out = table('Size', [height(db_in), width(db_in)], ...
                   'VariableTypes', repmat({'string'}, 1, width(db_in)), ...
                   'VariableNames', columnNames);
    
    for i = 1:width(db_in)
        db_out{:, i} = string(db_in{:, i});
    end
end
