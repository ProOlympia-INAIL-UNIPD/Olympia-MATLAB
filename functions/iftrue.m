function out = iftrue(in)
if ischar(in)
    if strcmp(in, 'True')
        out = 1;
    else
        out = 0;

    end
else
    out = in;
end