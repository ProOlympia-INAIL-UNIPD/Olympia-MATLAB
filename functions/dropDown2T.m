function T = dropDown2T(dval)
    T=eye(4);
    switch dval
        case 'Y-UP, Z-LAT'
        case 'Z-UP, X-LAT'
            T(:,[3 1 2 4])=T(:,[2 3 1 4]);
            
        case 'Z-UP, Y-LAT'
            %T(:,[3 1])=T(:,[2 3]);
    end
end