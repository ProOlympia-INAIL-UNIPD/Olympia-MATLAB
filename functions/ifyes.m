function out = ifyes(in)
if ischar(in)
    if strcmp(in, 'Yes')
        out = 1;
    else
        out = 0;

    end
else
    out = in;
end