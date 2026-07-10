function out = extractParameters(s, prefix)

if nargin < 2
    prefix = "";
end

out = struct;
fields = fieldnames(s);

for i = 1:length(fields)

    fname = fields{i};
    val   = s.(fname);

    % costruzione nome senza duplicati
    if strcmp(prefix, "")
        newname = fname;
    else
        
        parts = split(prefix,"_");
        last  = parts(end);

        if strcmp(last, fname)
            newname = prefix;
        else
            newname = prefix + "_" + fname;
        end
    end

    if isstruct(val)

        sub = extractParameters(val, newname);
        subfields = fieldnames(sub);

        for j = 1:length(subfields)
            out.(subfields{j}) = sub.(subfields{j});
        end

    else

        out.(newname) = val;

    end

end
end