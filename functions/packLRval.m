function arr = packLRval(nameL, nameR, context)
    if nargin < 3
        context = ["left","right"];
    end
    arr = [ struct("value", nameL, "context", context(1)), ...
            struct("value", nameR, "context", context(2)) ];
end