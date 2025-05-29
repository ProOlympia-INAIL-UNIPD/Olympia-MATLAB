function vec = zero2nan(vec)
%ZERO2NAN Summary of this function goes here
%   Detailed explanation goes here
vec(vec==0)=nan;
end

